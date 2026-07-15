#include "cch.h"
#include <algorithm>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

using namespace std;
using namespace RcppParallel;

namespace {

const double INF = numeric_limits<double>::max();

// This file implements the CCH (Customizable Contraction Hierarchy) backend.
//
// A CCH is split into three separate ideas:
// 1. Preparation: build a shortcut topology from the road graph alone.
// 2. Customization: put the current edge costs on that topology.
// 3. Query/AON: use the customized topology to find shortest paths and load
//    traffic flows back onto original road edges.
//
// This split is what makes CCH useful for congestion assignment. The road
// topology usually stays fixed while travel times change at every model
// iteration, so we pay the preparation cost once and only repeat customization
// and queries.

long long edge_key(int a, int b){
  // CCH preparation treats input edges as undirected. This key lets us find the
  // prepared CCH arc for either input direction of the same road link.
  if (a > b) std::swap(a, b);
  return (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);
}

IVec identity_order(int nb){
  IVec order(nb);
  for (int i = 0; i < nb; i++) order[i] = i;
  return order;
}

IVec min_degree_order(IVec &gfrom, IVec &gto, int nb){
  // Conservative built-in fallback. A lazy heap avoids scanning every node at
  // every elimination step. Large road networks should still pass a stronger
  // nested-dissection order because order quality dominates CCH size.
  vector<set<int> > graph(nb);
  for (int i = 0; i < gfrom.size(); i++){
    if (gfrom[i] == gto[i]) continue;
    graph[gfrom[i]].insert(gto[i]);
    graph[gto[i]].insert(gfrom[i]);
  }

  IVec order;
  order.reserve(nb);
  IVec removed(nb, 0);
  using DegreeNode = pair<size_t, int>;
  priority_queue<DegreeNode, vector<DegreeNode>, greater<DegreeNode>> queue;
  for (int node = 0; node < nb; node++) queue.push(make_pair(graph[node].size(), node));

  for (int step = 0; step < nb; step++){
    while (!queue.empty()){
      int node = queue.top().second;
      size_t degree = queue.top().first;
      if (removed[node] == 0 && degree == graph[node].size()) break;
      queue.pop();
    }
    if (queue.empty()) break;
    int best = queue.top().second;
    queue.pop();

    IVec neighbors;
    for (set<int>::iterator it = graph[best].begin(); it != graph[best].end(); ++it){
      if (removed[*it] == 0) neighbors.push_back(*it);
    }

    for (int i = 0; i < neighbors.size(); i++){
      for (int j = i + 1; j < neighbors.size(); j++){
        bool inserted = graph[neighbors[i]].insert(neighbors[j]).second;
        graph[neighbors[j]].insert(neighbors[i]);
        if (inserted){
          queue.push(make_pair(graph[neighbors[i]].size(), neighbors[i]));
          queue.push(make_pair(graph[neighbors[j]].size(), neighbors[j]));
        }
      }
    }

    for (int neighbor : neighbors){
      graph[neighbor].erase(best);
      queue.push(make_pair(graph[neighbor].size(), neighbor));
    }
    graph[best].clear();
    removed[best] = 1;
    order.push_back(best);
  }

  return order;
}

IVec rank_from_order(IVec order, int nb){
  if (order.size() == 0) order = identity_order(nb);
  if (order.size() != nb) throw invalid_argument("order length must be equal to the number of nodes");

  IVec rank(nb, -1);
  for (int i = 0; i < nb; i++){
    int node = order[i];
    if (node < 0 || node >= nb) throw invalid_argument("order contains an invalid node id");
    if (rank[node] != -1) throw invalid_argument("order contains duplicated node ids");
    rank[node] = i;
  }
  return rank;
}

void touch_node(int node, IVec &is_touched, IVec &touched){
  if (is_touched[node] == 0){
    is_touched[node] = 1;
    touched.push_back(node);
  }
}

void sort_unique(IVec &values){
  sort(values.begin(), values.end());
  values.erase(unique(values.begin(), values.end()), values.end());
}

IVec merge_unique_tail(IVec &base, IVec &added){
  IVec result;
  result.reserve(base.size() + added.size());
  merge(base.begin(), base.end(), added.begin(), added.end(), back_inserter(result));
  result.erase(unique(result.begin(), result.end()), result.end());
  return result;
}

int find_rank_arc(int tail_rank,
                  int head_rank,
                  IVec &rank_first_out,
                  IVec &rank_adj_head,
                  IVec &rank_adj_arc,
                  int &cursor){
  // During customization we scan many triangles of the form
  // lower_rank -> mid_rank -> upper_rank. Because adjacency is sorted, the
  // caller can reuse cursor instead of searching from the beginning each time.
  int end = rank_first_out[tail_rank + 1];
  while (cursor < end && rank_adj_head[cursor] < head_rank) cursor++;
  if (cursor < end && rank_adj_head[cursor] == head_rank) return rank_adj_arc[cursor];
  return -1;
}

IVec elimination_tree_from_up_graph(int nb, IVec &first_out, IVec &adj_head){
  // In the upward CCH graph, the first upward neighbor is the next ancestor in
  // the elimination tree. This is the structure used by RoutingKit's fast
  // customized queries.
  IVec parent(nb, -1);
  for (int node = 0; node < nb; node++){
    if (first_out[node] != first_out[node + 1]) parent[node] = adj_head[first_out[node]];
  }
  return parent;
}

template<class F>
void for_ancestors(IVec &parent, int node, const F &f){
  while (node != -1){
    if (!f(node)) return;
    node = parent[node];
  }
}

void touch_distance_node(int node, DVec &dist, IVec &pred_node, IVec &pred_arc, IVec &pred_dir,
                         IVec &is_touched, IVec &touched){
  if (is_touched[node] == 0){
    is_touched[node] = 1;
    touched.push_back(node);
  }
}

void relax_outgoing_elimination(int node,
                                IVec &first_out,
                                IVec &adj_head,
                                IVec &adj_arc,
                                DVec &weight,
                                int arc_dir,
                                DVec &dist,
                                IVec &pred_node,
                                IVec &pred_arc,
                                IVec &pred_dir,
                                IVec &is_touched,
                                IVec &touched){
  // Push distance labels upward away from the repeated endpoint. The predecessor
  // arrays remember enough information to later unpack shortcuts into original
  // road edges.
  if (dist[node] == INF) return;
  for (int i = first_out[node]; i < first_out[node + 1]; i++){
    int arc = adj_arc[i];
    double edge_weight = weight[arc];
    if (edge_weight == INF) continue;

    int next = adj_head[i];
    double tentative = dist[node] + edge_weight;
    if (tentative < dist[next]){
      touch_distance_node(next, dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      dist[next] = tentative;
      pred_node[next] = node;
      pred_arc[next] = arc;
      pred_dir[next] = arc_dir;
    }
  }
}

void relax_incoming_elimination(int node,
                                IVec &first_out,
                                IVec &adj_head,
                                IVec &adj_arc,
                                DVec &weight,
                                int arc_dir,
                                DVec &dist,
                                IVec &pred_node,
                                IVec &pred_arc,
                                IVec &pred_dir,
                                IVec &is_touched,
                                IVec &touched){
  // Pull distance labels back down toward the OD endpoints pinned in the same
  // group. This is the second half of the RoutingKit-style elimination query.
  for (int i = first_out[node]; i < first_out[node + 1]; i++){
    int arc = adj_arc[i];
    double edge_weight = weight[arc];
    if (edge_weight == INF) continue;

    int next = adj_head[i];
    if (dist[next] == INF) continue;
    double tentative = dist[next] + edge_weight;
    if (tentative < dist[node]){
      touch_distance_node(node, dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      dist[node] = tentative;
      pred_node[node] = next;
      pred_arc[node] = arc;
      pred_dir[node] = arc_dir;
    }
  }
}

void clear_elimination_state(DVec &dist,
                             IVec &pred_node,
                             IVec &pred_arc,
                             IVec &pred_dir,
                             IVec &is_touched,
                             IVec &touched){
  // Queries touch only a subset of nodes. Clearing just those nodes avoids an
  // O(number_of_nodes) reset after every source/destination group.
  for (int i = 0; i < touched.size(); i++){
    int node = touched[i];
    dist[node] = INF;
    pred_node[node] = -1;
    pred_arc[node] = -1;
    pred_dir[node] = -1;
    is_touched[node] = 0;
  }
  touched.clear();
}

void unpack_metric_arc(int arc,
                       int direction,
                       IVec &forward_first_arc,
                       IVec &forward_first_dir,
                       IVec &forward_second_arc,
                       IVec &forward_second_dir,
                       IVec &forward_original,
                       IVec &backward_first_arc,
                       IVec &backward_first_dir,
                       IVec &backward_second_arc,
                       IVec &backward_second_dir,
                       IVec &backward_original,
                       IVec &path_edges){
  // A customized CCH arc can be either an original road edge or a shortcut made
  // of two smaller CCH arcs. Traffic assignment needs flows on original road
  // edges, so path recovery recursively expands every shortcut before loading
  // demand.
  if (arc < 0) return;

  if (direction == 1){
    if (forward_original[arc] >= 0){
      path_edges.push_back(forward_original[arc]);
      return;
    }
    unpack_metric_arc(forward_first_arc[arc], forward_first_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
    unpack_metric_arc(forward_second_arc[arc], forward_second_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
  } else {
    if (backward_original[arc] >= 0){
      path_edges.push_back(backward_original[arc]);
      return;
    }
    unpack_metric_arc(backward_first_arc[arc], backward_first_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
    unpack_metric_arc(backward_second_arc[arc], backward_second_dir[arc],
                      forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                      backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                      path_edges);
  }
}

struct EndpointGroups {
  bool group_sources;
  IVec roots;
  vector<IVec> pairs;
};

EndpointGroups group_pairs_by_repeated_endpoint(const IVec &dep, const IVec &arr){
  unordered_set<int> sources(dep.begin(), dep.end());
  unordered_set<int> targets(arr.begin(), arr.end());

  EndpointGroups result;
  result.group_sources = sources.size() <= targets.size();

  const IVec &endpoint = result.group_sources ? dep : arr;
  unordered_map<int, int> root_to_group;
  size_t group_count = result.group_sources ? sources.size() : targets.size();
  root_to_group.reserve(group_count);
  result.roots.reserve(group_count);
  result.pairs.reserve(group_count);

  for (int pair = 0; pair < static_cast<int>(endpoint.size()); pair++){
    int root = endpoint[pair];
    auto inserted = root_to_group.emplace(root, static_cast<int>(result.pairs.size()));
    if (inserted.second){
      result.roots.push_back(root);
      result.pairs.push_back(IVec());
    }
    result.pairs[inserted.first->second].push_back(pair);
  }

  return result;
}

struct CCHDistancePairWorker : public Worker {
  // One-to-one distance queries on a customized CCH metric. This is useful for
  // distance APIs and correctness checks, but it does not recover paths and is
  // not the traffic-assignment path.
  int nb;
  IVec &first_out;
  IVec &adj_head;
  IVec &adj_arc;
  DVec &forward;
  DVec &backward;
  IVec dep;
  IVec arr;
  DVec &result;

  CCHDistancePairWorker(int nb,
                        IVec &first_out,
                        IVec &adj_head,
                        IVec &adj_arc,
                        DVec &forward,
                        DVec &backward,
                        IVec dep,
                        IVec arr,
                        DVec &result)
    : nb(nb), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      forward(forward), backward(backward), dep(dep), arr(arr), result(result) {}

  void operator()(std::size_t begin, std::size_t end){
    DVec dist_forward(nb, INF);
    DVec dist_backward(nb, INF);
    IVec seen_forward(nb, 0);
    IVec seen_backward(nb, 0);
    IVec is_touched(nb, 0);
    IVec touched;
    touched.reserve(1024);

    for (std::size_t k = begin; k != end; k++){
      int source = dep[k];
      int target = arr[k];

      dist_forward[source] = 0.0;
      dist_backward[target] = 0.0;
      seen_forward[source] = 1;
      seen_backward[target] = 1;
      touch_node(source, is_touched, touched);
      touch_node(target, is_touched, touched);

      PQ forward_queue;
      PQ backward_queue;
      forward_queue.push(make_pair(source, 0.0));
      backward_queue.push(make_pair(target, 0.0));

      double best = source == target ? 0.0 : INF;

      while (!forward_queue.empty() || !backward_queue.empty()){
        double forward_min = forward_queue.empty() ? INF : forward_queue.top().second;
        double backward_min = backward_queue.empty() ? INF : backward_queue.top().second;
        if (forward_min >= best && backward_min >= best) break;

        if (!forward_queue.empty() && forward_min <= backward_min){
          int node = forward_queue.top().first;
          double dist = forward_queue.top().second;
          forward_queue.pop();

          if (dist > dist_forward[node]) continue;

          if (seen_backward[node] == 1){
            double candidate = dist_forward[node] + dist_backward[node];
            if (candidate < best) best = candidate;
          }

          for (int i = first_out[node]; i < first_out[node + 1]; i++){
            int arc = adj_arc[i];
            double edge_weight = forward[arc];
            if (edge_weight == INF) continue;

            int next = adj_head[i];
            double tentative = dist + edge_weight;
            if (tentative < dist_forward[next]){
              dist_forward[next] = tentative;
              seen_forward[next] = 1;
              touch_node(next, is_touched, touched);
              forward_queue.push(make_pair(next, tentative));

              if (seen_backward[next] == 1){
                double candidate = tentative + dist_backward[next];
                if (candidate < best) best = candidate;
              }
            }
          }
        } else if (!backward_queue.empty()){
          int node = backward_queue.top().first;
          double dist = backward_queue.top().second;
          backward_queue.pop();

          if (dist > dist_backward[node]) continue;

          if (seen_forward[node] == 1){
            double candidate = dist_forward[node] + dist_backward[node];
            if (candidate < best) best = candidate;
          }

          for (int i = first_out[node]; i < first_out[node + 1]; i++){
            int arc = adj_arc[i];
            double edge_weight = backward[arc];
            if (edge_weight == INF) continue;

            int next = adj_head[i];
            double tentative = dist + edge_weight;
            if (tentative < dist_backward[next]){
              dist_backward[next] = tentative;
              seen_backward[next] = 1;
              touch_node(next, is_touched, touched);
              backward_queue.push(make_pair(next, tentative));

              if (seen_forward[next] == 1){
                double candidate = dist_forward[next] + tentative;
                if (candidate < best) best = candidate;
              }
            }
          }
        }
      }

      result[k] = best;

      for (int i = 0; i < touched.size(); i++){
        int node = touched[i];
        dist_forward[node] = INF;
        dist_backward[node] = INF;
        seen_forward[node] = 0;
        seen_backward[node] = 0;
        is_touched[node] = 0;
      }
      touched.clear();
    }
  }
};

} // namespace

CCHPrepared build_cch(IVec &gfrom, IVec &gto, int nb, IVec order){
  // Preparation phase.
  //
  // Input:
  // - original directed road edges, using integer node ids;
  // - an optional node order.
  //
  // Output:
  // - a directed upward shortcut graph;
  // - the elimination tree used by fast grouped queries;
  // - a map from each original road edge to the CCH arc where its cost enters.
  //
  // No travel-time/congestion weight is stored here. That is intentional: the
  // result can be saved and reused as long as the road topology is unchanged.
  if (order.size() == 0) order = min_degree_order(gfrom, gto, nb);

  CCHPrepared result;
  result.nbnode = nb;
  result.rank = rank_from_order(order, nb);

  // Build the upward graph induced by the order. Each input edge is first
  // treated as undirected; direction is restored during customization.
  vector<IVec> upward(nb);
  for (int i = 0; i < gfrom.size(); i++){
    int from_rank = result.rank[gfrom[i]];
    int to_rank = result.rank[gto[i]];
    if (from_rank == to_rank) continue;
    if (from_rank < to_rank){
      upward[from_rank].push_back(to_rank);
    } else {
      upward[to_rank].push_back(from_rank);
    }
  }

  for (int node_rank = 0; node_rank < nb; node_rank++) sort_unique(upward[node_rank]);

  vector<pair<int, int> > rank_arcs;
  rank_arcs.reserve(gfrom.size());

  for (int node_rank = 0; node_rank < nb; node_rank++){
    IVec &neighbors = upward[node_rank];
    if (neighbors.size() == 0) continue;

    for (int i = 0; i < neighbors.size(); i++) rank_arcs.push_back(make_pair(node_rank, neighbors[i]));

    // Eliminate this rank by making all higher-ranked neighbors adjacent. The
    // fill edges are what let later customization evaluate shortcut triangles.
    int lowest_upward_neighbor = neighbors[0];
    if (neighbors.size() > 1){
      IVec fill_edges(neighbors.begin() + 1, neighbors.end());
      upward[lowest_upward_neighbor] = merge_unique_tail(upward[lowest_upward_neighbor], fill_edges);
    }
  }

  sort(rank_arcs.begin(), rank_arcs.end());
  rank_arcs.erase(unique(rank_arcs.begin(), rank_arcs.end()), rank_arcs.end());

  result.tail.reserve(rank_arcs.size());
  result.head.reserve(rank_arcs.size());

  vector<IVec> cch_upward(nb);
  for (int arc = 0; arc < rank_arcs.size(); arc++){
    int lower_rank = rank_arcs[arc].first;
    int higher_rank = rank_arcs[arc].second;
    int lower_node = order[lower_rank];
    int higher_node = order[higher_rank];

    result.tail.push_back(lower_node);
    result.head.push_back(higher_node);
    cch_upward[lower_rank].push_back(higher_rank);
  }

  result.rank_first_out.assign(nb + 1, 0);
  for (int rank_id = 0; rank_id < nb; rank_id++) result.rank_first_out[rank_id + 1] = cch_upward[rank_id].size();
  for (int rank_id = 1; rank_id <= nb; rank_id++) result.rank_first_out[rank_id] += result.rank_first_out[rank_id - 1];

  result.rank_adj_head.resize(result.tail.size());
  result.rank_adj_arc.resize(result.tail.size());
  IVec rank_cursor = result.rank_first_out;
  for (int arc = 0; arc < rank_arcs.size(); arc++){
    int lower_rank = rank_arcs[arc].first;
    int higher_rank = rank_arcs[arc].second;
    int pos = rank_cursor[lower_rank]++;
    result.rank_adj_head[pos] = higher_rank;
    result.rank_adj_arc[pos] = arc;
  }

  result.first_out.assign(nb + 1, 0);
  for (int i = 0; i < result.tail.size(); i++) result.first_out[result.tail[i] + 1]++;
  for (int i = 1; i <= nb; i++) result.first_out[i] += result.first_out[i - 1];

  result.adj_head.resize(result.tail.size());
  result.adj_arc.resize(result.tail.size());
  IVec cursor = result.first_out;
  for (int arc = 0; arc < result.tail.size(); arc++){
    int pos = cursor[result.tail[arc]]++;
    result.adj_head[pos] = result.head[arc];
    result.adj_arc[pos] = arc;
  }
  result.elimination_tree_parent = elimination_tree_from_up_graph(nb, result.first_out, result.adj_head);

  // Remember where each original edge enters the CCH. Congestion assignment
  // changes edge costs many times, so customization must not search for this
  // mapping at every iteration.
  unordered_map<long long, int> edge_to_arc;
  edge_to_arc.reserve(result.tail.size() * 2);
  for (int arc = 0; arc < result.tail.size(); arc++){
    edge_to_arc[edge_key(result.tail[arc], result.head[arc])] = arc;
  }

  result.input_arc.assign(gfrom.size(), -1);
  result.input_forward.assign(gfrom.size(), 0);
  for (int i = 0; i < gfrom.size(); i++){
    int from = gfrom[i];
    int to = gto[i];
    if (from == to) continue;

    int lower = result.rank[from] < result.rank[to] ? from : to;
    int higher = result.rank[from] < result.rank[to] ? to : from;
    unordered_map<long long, int>::iterator found = edge_to_arc.find(edge_key(lower, higher));
    if (found == edge_to_arc.end()) continue;

    result.input_arc[i] = found->second;
    result.input_forward[i] = (from == lower && to == higher) ? 1 : 0;
  }

  return result;
}

void remap_cch_input_arcs(IVec &gfrom, IVec &gto, CCHPrepared &prepared){
  remap_cch_input_arcs(gfrom, gto, prepared, prepared.input_arc, prepared.input_forward);
}

void remap_cch_input_arcs(IVec &gfrom,
                          IVec &gto,
                          CCHPrepared &prepared,
                          IVec &input_arc,
                          IVec &input_forward){
  // Prepared CCH objects can be reused by assign_traffic(). The assignment
  // graph has the same topology but its current adjacency edge order may differ,
  // so only this cheap input-edge map is rebuilt.
  unordered_map<long long, int> edge_to_arc;
  edge_to_arc.reserve(prepared.tail.size() * 2);
  for (int arc = 0; arc < prepared.tail.size(); arc++){
    edge_to_arc[edge_key(prepared.tail[arc], prepared.head[arc])] = arc;
  }

  input_arc.assign(gfrom.size(), -1);
  input_forward.assign(gfrom.size(), 0);
  for (int i = 0; i < gfrom.size(); i++){
    int from = gfrom[i];
    int to = gto[i];
    if (from == to) continue;

    int lower = prepared.rank[from] < prepared.rank[to] ? from : to;
    int higher = prepared.rank[from] < prepared.rank[to] ? to : from;
    unordered_map<long long, int>::iterator found = edge_to_arc.find(edge_key(lower, higher));
    if (found == edge_to_arc.end()) continue;

    input_arc[i] = found->second;
    input_forward[i] = (from == lower && to == higher) ? 1 : 0;
  }
}

void customize_cch(IVec &gfrom,
                   IVec &gto,
                   DVec &gw,
                   int nb,
                   IVec &rank,
                   IVec &tail,
                   IVec &head,
                   IVec &rank_first_out,
                   IVec &rank_adj_head,
                   IVec &rank_adj_arc,
                   IVec &input_arc,
                   IVec &input_forward,
                   DVec &forward,
                   DVec &backward,
                   IVec &forward_first_arc,
                   IVec &forward_first_dir,
                   IVec &forward_second_arc,
                   IVec &forward_second_dir,
                   IVec &forward_original,
                   IVec &backward_first_arc,
                   IVec &backward_first_dir,
                   IVec &backward_second_arc,
                   IVec &backward_second_dir,
                   IVec &backward_original){
  // Customization phase.
  //
  // The caller supplies one current cost per original edge. This function
  // updates the prepared shortcut graph so every CCH arc stores the best known
  // cost in each travel direction:
  // - forward: lower-ranked node -> higher-ranked node;
  // - backward: higher-ranked node -> lower-ranked node.
  //
  // For original arcs, *_original points to the input edge. For shortcuts,
  // *_first_arc and *_second_arc describe how to unpack the shortcut later.
  // CCH customization has two jobs:
  // 1. Load current edge weights onto the prepared CCH arcs.
  // 2. Relax all upward triangles so each shortcut stores the cheapest unpackable
  //    path for the current congestion state.
  forward.assign(tail.size(), INF);
  backward.assign(tail.size(), INF);
  forward_first_arc.assign(tail.size(), -1);
  forward_first_dir.assign(tail.size(), -1);
  forward_second_arc.assign(tail.size(), -1);
  forward_second_dir.assign(tail.size(), -1);
  forward_original.assign(tail.size(), -1);
  backward_first_arc.assign(tail.size(), -1);
  backward_first_dir.assign(tail.size(), -1);
  backward_second_arc.assign(tail.size(), -1);
  backward_second_dir.assign(tail.size(), -1);
  backward_original.assign(tail.size(), -1);

  for (int i = 0; i < gfrom.size(); i++){
    int arc = input_arc[i];
    if (arc < 0) continue;

    if (input_forward[i] == 1){
      if (gw[i] < forward[arc]){
        forward[arc] = gw[i];
        forward_original[arc] = i;
        forward_first_arc[arc] = -1;
        forward_second_arc[arc] = -1;
      }
    } else {
      if (gw[i] < backward[arc]){
        backward[arc] = gw[i];
        backward_original[arc] = i;
        backward_first_arc[arc] = -1;
        backward_second_arc[arc] = -1;
      }
    }
  }

  for (int lower_rank = 0; lower_rank < nb; lower_rank++){
    int begin = rank_first_out[lower_rank];
    int end = rank_first_out[lower_rank + 1];

    for (int i = begin; i < end; i++){
      int mid_rank = rank_adj_head[i];
      int left = rank_adj_arc[i];
      int mid_cursor = rank_first_out[mid_rank];

      for (int j = i + 1; j < end; j++){
        int upper_rank = rank_adj_head[j];
        int right = rank_adj_arc[j];
        int upper = find_rank_arc(mid_rank, upper_rank, rank_first_out, rank_adj_head, rank_adj_arc, mid_cursor);
        if (upper < 0) continue;

        // left and right form lower -> mid -> upper in rank space. Depending on
        // the original travel direction, this can improve either metric
        // direction stored on the upper shortcut.
        if (backward[left] != INF && forward[right] != INF){
          double candidate = backward[left] + forward[right];
          if (candidate < forward[upper]){
            forward[upper] = candidate;
            forward_original[upper] = -1;
            forward_first_arc[upper] = left;
            forward_first_dir[upper] = 0;
            forward_second_arc[upper] = right;
            forward_second_dir[upper] = 1;
          }
        }

        if (backward[right] != INF && forward[left] != INF){
          double candidate = backward[right] + forward[left];
          if (candidate < backward[upper]){
            backward[upper] = candidate;
            backward_original[upper] = -1;
            backward_first_arc[upper] = right;
            backward_first_dir[upper] = 0;
            backward_second_arc[upper] = left;
            backward_second_dir[upper] = 1;
          }
        }
      }
    }
  }
}

DVec distance_pair_cch(int nb,
                       IVec &first_out,
                       IVec &adj_head,
                       IVec &adj_arc,
                       DVec &forward,
                       DVec &backward,
                       IVec dep,
                       IVec arr){
  // Distance-only pair query. It computes shortest path costs but intentionally
  // does not recover original road edges, so it cannot be used for AON flow
  // loading.
  DVec result(dep.size(), INF);
  CCHDistancePairWorker worker(nb, first_out, adj_head, adj_arc, forward, backward, dep, arr, result);
  if (dep.size() <= 1000){
    worker(0, dep.size());
  } else {
    parallelFor(0, dep.size(), worker);
  }

  return result;
}

struct CCHDistanceMatrixWorker : public Worker {
  int nb;
  IVec &first_out;
  IVec &adj_head;
  IVec &adj_arc;
  DVec &outgoing_weight;
  DVec &incoming_weight;
  IVec &elimination_tree_parent;
  IVec &roots;
  IVec &pinned;
  IVec &pinned_stack;
  bool roots_are_sources;
  RMatrix<double> result;

  CCHDistanceMatrixWorker(int nb,
                          IVec &first_out,
                          IVec &adj_head,
                          IVec &adj_arc,
                          DVec &outgoing_weight,
                          DVec &incoming_weight,
                          IVec &elimination_tree_parent,
                          IVec &roots,
                          IVec &pinned,
                          IVec &pinned_stack,
                          bool roots_are_sources,
                          Rcpp::NumericMatrix result)
    : nb(nb), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      outgoing_weight(outgoing_weight), incoming_weight(incoming_weight),
      elimination_tree_parent(elimination_tree_parent), roots(roots), pinned(pinned),
      pinned_stack(pinned_stack), roots_are_sources(roots_are_sources), result(result) {}

  void relax_outgoing(int node, DVec &dist, IVec &is_touched, IVec &touched){
    if (dist[node] == INF) return;
    for (int i = first_out[node]; i < first_out[node + 1]; i++){
      int arc = adj_arc[i];
      double weight = outgoing_weight[arc];
      if (weight == INF) continue;
      int next = adj_head[i];
      double candidate = dist[node] + weight;
      if (candidate < dist[next]){
        touch_node(next, is_touched, touched);
        dist[next] = candidate;
      }
    }
  }

  void relax_incoming(int node, DVec &dist, IVec &is_touched, IVec &touched){
    for (int i = first_out[node]; i < first_out[node + 1]; i++){
      int arc = adj_arc[i];
      double weight = incoming_weight[arc];
      int next = adj_head[i];
      if (weight == INF || dist[next] == INF) continue;
      double candidate = dist[next] + weight;
      if (candidate < dist[node]){
        touch_node(node, is_touched, touched);
        dist[node] = candidate;
      }
    }
  }

  void operator()(std::size_t begin, std::size_t end){
    DVec dist(nb, INF);
    IVec is_touched(nb, 0);
    IVec touched;

    for (std::size_t root_index = begin; root_index != end; root_index++){
      int root = roots[root_index];
      touch_node(root, is_touched, touched);
      dist[root] = 0.0;

      for_ancestors(elimination_tree_parent, root, [&](int node){
        relax_outgoing(node, dist, is_touched, touched);
        return true;
      });

      for (int node : pinned_stack) relax_incoming(node, dist, is_touched, touched);

      for (int pinned_index = 0; pinned_index < static_cast<int>(pinned.size()); pinned_index++){
        if (roots_are_sources){
          result(root_index, pinned_index) = dist[pinned[pinned_index]];
        } else {
          result(pinned_index, root_index) = dist[pinned[pinned_index]];
        }
      }

      for (int node : touched){
        dist[node] = INF;
        is_touched[node] = 0;
      }
      touched.clear();
    }
  }
};

Rcpp::NumericMatrix distance_matrix_cch(int nb,
                                        IVec &rank,
                                        IVec &first_out,
                                        IVec &adj_head,
                                        IVec &adj_arc,
                                        IVec &elimination_tree_parent,
                                        DVec &forward,
                                        DVec &backward,
                                        IVec dep,
                                        IVec arr){
  Rcpp::NumericMatrix result(dep.size(), arr.size());
  if (dep.empty() || arr.empty()) return result;

  bool roots_are_sources = dep.size() <= arr.size();
  IVec &roots = roots_are_sources ? dep : arr;
  IVec &pinned = roots_are_sources ? arr : dep;
  DVec &outgoing_weight = roots_are_sources ? forward : backward;
  DVec &incoming_weight = roots_are_sources ? backward : forward;

  IVec is_pinned_ancestor(nb, 0);
  IVec pinned_stack;
  for (int endpoint : pinned){
    for_ancestors(elimination_tree_parent, endpoint, [&](int node){
      if (is_pinned_ancestor[node] == 0){
        is_pinned_ancestor[node] = 1;
        pinned_stack.push_back(node);
      }
      return true;
    });
  }
  sort(pinned_stack.begin(), pinned_stack.end(), [&](int a, int b){
    return rank[a] > rank[b];
  });

  CCHDistanceMatrixWorker worker(nb, first_out, adj_head, adj_arc,
                                 outgoing_weight, incoming_weight,
                                 elimination_tree_parent, roots, pinned,
                                 pinned_stack, roots_are_sources, result);
  if (roots.size() <= 1 || nb <= 1000){
    worker(0, roots.size());
  } else {
    parallelFor(0, roots.size(), worker);
  }

  return result;
}

struct CCHPathValuesPairWorker : public Worker {
  // Pair query that routes on the customized CCH cost and accumulates one or
  // more original-edge value columns along the selected shortest path.
  //
  // This is the Mobility-shaped operation: keep one routing cost, then collect
  // distance, real time, leg distances, or any other edge attribute in one pass.
  IVec &gfrom;
  int nb;
  int n_values;
  IVec &rank;
  IVec &first_out;
  IVec &adj_head;
  IVec &adj_arc;
  IVec &elimination_tree_parent;
  DVec &forward;
  DVec &backward;
  IVec &forward_first_arc;
  IVec &forward_first_dir;
  IVec &forward_second_arc;
  IVec &forward_second_dir;
  IVec &forward_original;
  IVec &backward_first_arc;
  IVec &backward_first_dir;
  IVec &backward_second_arc;
  IVec &backward_second_dir;
  IVec &backward_original;
  IVec &dep;
  IVec &arr;
  DVec &values;
  IVec &group_roots;
  vector<IVec> &groups;
  bool group_sources;
  DVec &result;

  CCHPathValuesPairWorker(IVec &gfrom,
                          int nb,
                          int n_values,
                          IVec &rank,
                          IVec &first_out,
                          IVec &adj_head,
                          IVec &adj_arc,
                          IVec &elimination_tree_parent,
                          DVec &forward,
                          DVec &backward,
                          IVec &forward_first_arc,
                          IVec &forward_first_dir,
                          IVec &forward_second_arc,
                          IVec &forward_second_dir,
                          IVec &forward_original,
                          IVec &backward_first_arc,
                          IVec &backward_first_dir,
                          IVec &backward_second_arc,
                          IVec &backward_second_dir,
                          IVec &backward_original,
                          IVec &dep,
                          IVec &arr,
                          DVec &values,
                          IVec &group_roots,
                          vector<IVec> &groups,
                          bool group_sources,
                          DVec &result)
    : gfrom(gfrom), nb(nb), n_values(n_values), rank(rank), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      elimination_tree_parent(elimination_tree_parent), forward(forward), backward(backward),
      forward_first_arc(forward_first_arc), forward_first_dir(forward_first_dir),
      forward_second_arc(forward_second_arc), forward_second_dir(forward_second_dir),
      forward_original(forward_original),
      backward_first_arc(backward_first_arc), backward_first_dir(backward_first_dir),
      backward_second_arc(backward_second_arc), backward_second_dir(backward_second_dir),
      backward_original(backward_original), dep(dep), arr(arr), values(values),
      group_roots(group_roots), groups(groups), group_sources(group_sources), result(result) {}

  CCHPathValuesPairWorker(CCHPathValuesPairWorker &other, Split)
    : gfrom(other.gfrom), nb(other.nb), n_values(other.n_values), rank(other.rank),
      first_out(other.first_out), adj_head(other.adj_head), adj_arc(other.adj_arc),
      elimination_tree_parent(other.elimination_tree_parent),
      forward(other.forward), backward(other.backward),
      forward_first_arc(other.forward_first_arc), forward_first_dir(other.forward_first_dir),
      forward_second_arc(other.forward_second_arc), forward_second_dir(other.forward_second_dir),
      forward_original(other.forward_original),
      backward_first_arc(other.backward_first_arc), backward_first_dir(other.backward_first_dir),
      backward_second_arc(other.backward_second_arc), backward_second_dir(other.backward_second_dir),
      backward_original(other.backward_original), dep(other.dep), arr(other.arr), values(other.values),
      group_roots(other.group_roots), groups(other.groups), group_sources(other.group_sources), result(other.result) {}

  void add_recovered_values(int pair_index,
                            int start,
                            int root,
                            bool reverse_path,
                            double cost,
                            IVec &pred_node,
                            IVec &pred_arc,
                            IVec &pred_dir,
                            IVec &path_arcs,
                            IVec &path_dirs,
                            IVec &path_edges){
    int n_pairs = dep.size();
    result[pair_index] = cost;
    for (int col = 0; col < n_values; col++) result[(col + 1) * n_pairs + pair_index] = 0.0;

    if (start == root) return;

    int node = start;
    while (node != root && pred_node[node] != -1){
      path_arcs.push_back(pred_arc[node]);
      path_dirs.push_back(pred_dir[node]);
      node = pred_node[node];
    }
    if (node != root){
      result[pair_index] = INF;
      for (int col = 0; col < n_values; col++) result[(col + 1) * n_pairs + pair_index] = INF;
      return;
    }
    if (reverse_path){
      reverse(path_arcs.begin(), path_arcs.end());
      reverse(path_dirs.begin(), path_dirs.end());
    }

    for (int i = 0; i < path_arcs.size(); i++){
      unpack_metric_arc(path_arcs[i], path_dirs[i],
                        forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                        backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                        path_edges);
    }
    for (int i = 0; i < path_edges.size(); i++){
      int edge = path_edges[i];
      if (edge < 0) continue;
      for (int col = 0; col < n_values; col++){
        result[(col + 1) * n_pairs + pair_index] += values[col * gfrom.size() + edge];
      }
    }
  }

  void operator()(std::size_t begin, std::size_t end){
    DVec dist(nb, INF);
    IVec pred_node(nb, -1);
    IVec pred_arc(nb, -1);
    IVec pred_dir(nb, -1);
    IVec is_touched(nb, 0);
    IVec touched;
    IVec is_pinned_ancestor(nb, 0);
    IVec pinned_stack;
    IVec pinned_touched;
    IVec path_arcs;
    IVec path_dirs;
    IVec path_edges;

    for (std::size_t group_id = begin; group_id != end; group_id++){
      if (groups[group_id].size() == 0) continue;

      int root = group_roots[group_id];
      touch_distance_node(root, dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      dist[root] = 0.0;

      if (group_sources){
        for_ancestors(elimination_tree_parent, root, [&](int node){
          relax_outgoing_elimination(node, first_out, adj_head, adj_arc, forward, 1,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
          return true;
        });
      } else {
        for_ancestors(elimination_tree_parent, root, [&](int node){
          relax_outgoing_elimination(node, first_out, adj_head, adj_arc, backward, 0,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
          return true;
        });
      }

      for (int pos = 0; pos < groups[group_id].size(); pos++){
        int k = groups[group_id][pos];
        int pinned = group_sources ? arr[k] : dep[k];
        for_ancestors(elimination_tree_parent, pinned, [&](int node){
          if (is_pinned_ancestor[node] == 0){
            is_pinned_ancestor[node] = 1;
            pinned_touched.push_back(node);
            pinned_stack.push_back(node);
          }
          return true;
        });
      }

      sort(pinned_stack.begin(), pinned_stack.end(), [&](int a, int b){
        return rank[a] > rank[b];
      });
      for (int stack_pos = 0; stack_pos < pinned_stack.size(); stack_pos++){
        int node = pinned_stack[stack_pos];
        if (group_sources){
          relax_incoming_elimination(node, first_out, adj_head, adj_arc, backward, 0,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
        } else {
          relax_incoming_elimination(node, first_out, adj_head, adj_arc, forward, 1,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
        }
      }

      for (int pos = 0; pos < groups[group_id].size(); pos++){
        int k = groups[group_id][pos];
        int start = group_sources ? arr[k] : dep[k];
        if (dist[start] == INF) continue;

        path_arcs.clear();
        path_dirs.clear();
        path_edges.clear();
        add_recovered_values(k, start, root, group_sources, dist[start],
                             pred_node, pred_arc, pred_dir, path_arcs, path_dirs, path_edges);
      }

      clear_elimination_state(dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      for (int i = 0; i < pinned_touched.size(); i++) is_pinned_ancestor[pinned_touched[i]] = 0;
      pinned_touched.clear();
      pinned_stack.clear();
    }
  }
};

Rcpp::NumericMatrix path_values_pair_cch(IVec &gfrom,
                                         int nb,
                                         IVec &rank,
                                         IVec &first_out,
                                         IVec &adj_head,
                                         IVec &adj_arc,
                                         IVec &elimination_tree_parent,
                                         DVec &forward,
                                         DVec &backward,
                                         IVec &forward_first_arc,
                                         IVec &forward_first_dir,
                                         IVec &forward_second_arc,
                                         IVec &forward_second_dir,
                                         IVec &forward_original,
                                         IVec &backward_first_arc,
                                         IVec &backward_first_dir,
                                         IVec &backward_second_arc,
                                         IVec &backward_second_dir,
                                         IVec &backward_original,
                                         IVec dep,
                                         IVec arr,
                                         DVec &values,
                                         int n_values){
  EndpointGroups endpoint_groups = group_pairs_by_repeated_endpoint(dep, arr);

  DVec result((n_values + 1) * dep.size(), INF);
  CCHPathValuesPairWorker worker(gfrom, nb, n_values, rank, first_out, adj_head, adj_arc, elimination_tree_parent,
                                 forward, backward,
                                 forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                                 backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                                 dep, arr, values, endpoint_groups.roots, endpoint_groups.pairs,
                                 endpoint_groups.group_sources, result);
  if (dep.size() <= 1000 || nb <= 1000 || endpoint_groups.pairs.size() <= 1){
    worker(0, endpoint_groups.pairs.size());
  } else {
    parallelFor(0, endpoint_groups.pairs.size(), worker);
  }

  Rcpp::NumericMatrix out(dep.size(), n_values + 1);
  for (int col = 0; col < n_values + 1; col++){
    for (int row = 0; row < dep.size(); row++){
      out(row, col) = result[col * dep.size() + row];
    }
  }
  return out;
}

struct CCHEliminationGroupedAonWorker : public Worker {
  // RoutingKit-style grouped elimination-tree AON query.
  //
  // Instead of running one independent shortest-path search per OD pair, this
  // worker groups OD pairs by the side with fewer unique endpoints. A group then
  // shares one search from that repeated endpoint and only performs the smaller
  // pinned-endpoint work needed for the individual OD pairs.
  //
  // The recovered paths are unpacked to original graph edges immediately and
  // accumulated in flow. That is why this worker is the traffic-assignment CCH
  // path, while the distance query above is not.
  IVec &gfrom;
  int nb;
  IVec &rank;
  IVec &first_out;
  IVec &adj_head;
  IVec &adj_arc;
  IVec &elimination_tree_parent;
  DVec &forward;
  DVec &backward;
  IVec &forward_first_arc;
  IVec &forward_first_dir;
  IVec &forward_second_arc;
  IVec &forward_second_dir;
  IVec &forward_original;
  IVec &backward_first_arc;
  IVec &backward_first_dir;
  IVec &backward_second_arc;
  IVec &backward_second_dir;
  IVec &backward_original;
  IVec &dep;
  IVec &arr;
  DVec &demand;
  IVec &group_roots;
  vector<IVec> &groups;
  bool group_sources;
  DVec flow;

  CCHEliminationGroupedAonWorker(IVec &gfrom,
                                 int nb,
                                 IVec &rank,
                                 IVec &first_out,
                                 IVec &adj_head,
                                 IVec &adj_arc,
                                 IVec &elimination_tree_parent,
                                 DVec &forward,
                                 DVec &backward,
                                 IVec &forward_first_arc,
                                 IVec &forward_first_dir,
                                 IVec &forward_second_arc,
                                 IVec &forward_second_dir,
                                 IVec &forward_original,
                                 IVec &backward_first_arc,
                                 IVec &backward_first_dir,
                                 IVec &backward_second_arc,
                                 IVec &backward_second_dir,
                                 IVec &backward_original,
                                 IVec &dep,
                                 IVec &arr,
                                 DVec &demand,
                                 IVec &group_roots,
                                 vector<IVec> &groups,
                                 bool group_sources)
    : gfrom(gfrom), nb(nb), rank(rank), first_out(first_out), adj_head(adj_head), adj_arc(adj_arc),
      elimination_tree_parent(elimination_tree_parent), forward(forward), backward(backward),
      forward_first_arc(forward_first_arc), forward_first_dir(forward_first_dir),
      forward_second_arc(forward_second_arc), forward_second_dir(forward_second_dir),
      forward_original(forward_original),
      backward_first_arc(backward_first_arc), backward_first_dir(backward_first_dir),
      backward_second_arc(backward_second_arc), backward_second_dir(backward_second_dir),
      backward_original(backward_original), dep(dep), arr(arr), demand(demand),
      group_roots(group_roots), groups(groups), group_sources(group_sources), flow(gfrom.size(), 0.0) {}

  CCHEliminationGroupedAonWorker(CCHEliminationGroupedAonWorker &other, Split)
    : gfrom(other.gfrom), nb(other.nb), rank(other.rank),
      first_out(other.first_out), adj_head(other.adj_head), adj_arc(other.adj_arc),
      elimination_tree_parent(other.elimination_tree_parent),
      forward(other.forward), backward(other.backward),
      forward_first_arc(other.forward_first_arc), forward_first_dir(other.forward_first_dir),
      forward_second_arc(other.forward_second_arc), forward_second_dir(other.forward_second_dir),
      forward_original(other.forward_original),
      backward_first_arc(other.backward_first_arc), backward_first_dir(other.backward_first_dir),
      backward_second_arc(other.backward_second_arc), backward_second_dir(other.backward_second_dir),
      backward_original(other.backward_original), dep(other.dep), arr(other.arr), demand(other.demand),
      group_roots(other.group_roots), groups(other.groups), group_sources(other.group_sources), flow(other.gfrom.size(), 0.0) {}

  void join(CCHEliminationGroupedAonWorker &other){
    for (int i = 0; i < flow.size(); i++) flow[i] += other.flow[i];
  }

  void add_recovered_path(int start,
                          int root,
                          bool reverse_path,
                          double volume,
                          IVec &pred_node,
                          IVec &pred_arc,
                          IVec &pred_dir,
                          IVec &path_arcs,
                          IVec &path_dirs,
                          IVec &path_edges){
    // pred_* describes the path from start back to the group root in CCH arc
    // space. reverse_path is needed because source-group and destination-group
    // recovery walk the same predecessor tree in opposite logical directions.
    int node = start;
    while (node != root && pred_node[node] != -1){
      path_arcs.push_back(pred_arc[node]);
      path_dirs.push_back(pred_dir[node]);
      node = pred_node[node];
    }
    if (node != root) {
      path_arcs.clear();
      path_dirs.clear();
      return;
    }
    if (reverse_path){
      reverse(path_arcs.begin(), path_arcs.end());
      reverse(path_dirs.begin(), path_dirs.end());
    }

    for (int i = 0; i < path_arcs.size(); i++){
      unpack_metric_arc(path_arcs[i], path_dirs[i],
                        forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                        backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                        path_edges);
    }
    for (int i = 0; i < path_edges.size(); i++){
      if (path_edges[i] >= 0) flow[path_edges[i]] += volume;
    }
  }

  void operator()(std::size_t begin, std::size_t end){
    DVec dist(nb, INF);
    IVec pred_node(nb, -1);
    IVec pred_arc(nb, -1);
    IVec pred_dir(nb, -1);
    IVec is_touched(nb, 0);
    IVec touched;
    IVec is_pinned_ancestor(nb, 0);
    IVec pinned_stack;
    IVec pinned_touched;
    IVec path_arcs;
    IVec path_dirs;
    IVec path_edges;

    for (std::size_t group_id = begin; group_id != end; group_id++){
      if (groups[group_id].size() == 0) continue;

      int root = group_roots[group_id];
      touch_distance_node(root, dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      dist[root] = 0.0;

      // Shared phase: one upward search from the repeated endpoint. For source
      // groups this is a forward search from the origin; for destination groups
      // it is the symmetric backward search from the destination.
      if (group_sources){
        for_ancestors(elimination_tree_parent, root, [&](int node){
          relax_outgoing_elimination(node, first_out, adj_head, adj_arc, forward, 1,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
          return true;
        });
      } else {
        for_ancestors(elimination_tree_parent, root, [&](int node){
          relax_outgoing_elimination(node, first_out, adj_head, adj_arc, backward, 0,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
          return true;
        });
      }

      // Pinned phase: only visit ancestors of the OD endpoints in this group.
      // Processing high ranks first is essential. It makes the dynamic program
      // propagate distances from common ancestors back down toward each pinned
      // endpoint before we recover paths.
      for (int pos = 0; pos < groups[group_id].size(); pos++){
        int k = groups[group_id][pos];
        int pinned = group_sources ? arr[k] : dep[k];
        for_ancestors(elimination_tree_parent, pinned, [&](int node){
          if (is_pinned_ancestor[node] == 0){
            is_pinned_ancestor[node] = 1;
            pinned_touched.push_back(node);
            pinned_stack.push_back(node);
          }
          return true;
        });
      }

      sort(pinned_stack.begin(), pinned_stack.end(), [&](int a, int b){
        return rank[a] > rank[b];
      });
      for (int stack_pos = 0; stack_pos < pinned_stack.size(); stack_pos++){
        int node = pinned_stack[stack_pos];
        if (group_sources){
          relax_incoming_elimination(node, first_out, adj_head, adj_arc, backward, 0,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
        } else {
          relax_incoming_elimination(node, first_out, adj_head, adj_arc, forward, 1,
                                     dist, pred_node, pred_arc, pred_dir, is_touched, touched);
        }
      }

      // Recover and unpack one shortest path for every OD pair in the group.
      // Shortcut unpacking writes flow back to original input edges, which is
      // what traffic assignment needs for the next cost update.
      for (int pos = 0; pos < groups[group_id].size(); pos++){
        int k = groups[group_id][pos];
        int start = group_sources ? arr[k] : dep[k];
        if (start == root || dist[start] == INF) continue;

        path_arcs.clear();
        path_dirs.clear();
        path_edges.clear();
        add_recovered_path(start, root, group_sources, demand[k],
                           pred_node, pred_arc, pred_dir, path_arcs, path_dirs, path_edges);
      }

      clear_elimination_state(dist, pred_node, pred_arc, pred_dir, is_touched, touched);
      for (int i = 0; i < pinned_touched.size(); i++) is_pinned_ancestor[pinned_touched[i]] = 0;
      pinned_touched.clear();
      pinned_stack.clear();
    }
  }
};

DVec aon_flow_cch_elimination_grouped(IVec &gfrom,
                                      IVec &gto,
                                      DVec &gw,
                                      int nb,
                                      IVec &rank,
                                      IVec &first_out,
                                      IVec &adj_head,
                                      IVec &adj_arc,
                                      IVec &elimination_tree_parent,
                                      DVec &forward,
                                      DVec &backward,
                                      IVec &forward_first_arc,
                                      IVec &forward_first_dir,
                                      IVec &forward_second_arc,
                                      IVec &forward_second_dir,
                                      IVec &forward_original,
                                      IVec &backward_first_arc,
                                      IVec &backward_first_dir,
                                      IVec &backward_second_arc,
                                      IVec &backward_second_dir,
                                      IVec &backward_original,
                                      IVec &dep,
                                      IVec &arr,
                                      DVec &demand){
  // Public native AON entry point for the current CCH implementation. The R API
  // exposes this as aon_method = "cch". Older pairwise/grouped CCH query modes
  // were removed so assignment always uses this elimination-tree path.
  EndpointGroups endpoint_groups = group_pairs_by_repeated_endpoint(dep, arr);

  CCHEliminationGroupedAonWorker worker(gfrom, nb, rank, first_out, adj_head, adj_arc, elimination_tree_parent,
                                        forward, backward,
                                        forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                                        backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                                        dep, arr, demand, endpoint_groups.roots, endpoint_groups.pairs,
                                        endpoint_groups.group_sources);
  if (dep.size() <= 1000 || nb <= 1000 || endpoint_groups.pairs.size() <= 1){
    worker(0, endpoint_groups.pairs.size());
  } else {
    parallelReduce(0, endpoint_groups.pairs.size(), worker);
  }

  return worker.flow;
}

Rcpp::List aon_cch_elimination_grouped(IVec &gfrom,
                                       IVec &gto,
                                       DVec &gw,
                                       int nb,
                                       IVec &rank,
                                       IVec &first_out,
                                       IVec &adj_head,
                                       IVec &adj_arc,
                                       IVec &elimination_tree_parent,
                                       DVec &forward,
                                       DVec &backward,
                                       IVec &forward_first_arc,
                                       IVec &forward_first_dir,
                                       IVec &forward_second_arc,
                                       IVec &forward_second_dir,
                                       IVec &forward_original,
                                       IVec &backward_first_arc,
                                       IVec &backward_first_dir,
                                       IVec &backward_second_arc,
                                       IVec &backward_second_dir,
                                       IVec &backward_original,
                                       IVec dep,
                                       IVec arr,
                                       DVec demand){
  // Thin Rcpp-friendly wrapper around the native flow vector. It returns the
  // same four-column shape as the other cppRouting AON backends.
  DVec flow = aon_flow_cch_elimination_grouped(gfrom, gto, gw, nb, rank, first_out, adj_head, adj_arc,
                                               elimination_tree_parent, forward, backward,
                                               forward_first_arc, forward_first_dir, forward_second_arc, forward_second_dir, forward_original,
                                               backward_first_arc, backward_first_dir, backward_second_arc, backward_second_dir, backward_original,
                                               dep, arr, demand);

  Rcpp::List result(4);
  result[0] = gfrom;
  result[1] = gto;
  result[2] = gw;
  result[3] = flow;
  return result;
}
