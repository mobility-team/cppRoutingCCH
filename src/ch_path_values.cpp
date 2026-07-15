#include "cgraph.h"
#include "stall.h"
#include <Rcpp.h>
#include <RcppParallel.h>
#include <limits>
#include <algorithm>
#include <vector>

using namespace RcppParallel;
using namespace std;

using OriginalEdgeIndex = vector<vector<pair<int, int>>>;

struct CHPathValuesPairWorker : public Worker {
  CGraph* graph;
  const IVec& dep;
  const IVec& arr;
  const DVec& original_weight;
  const OriginalEdgeIndex& original_edges;
  RMatrix<double> values;
  RMatrix<double> result;

  CHPathValuesPairWorker(CGraph* graph,
                         const IVec& dep,
                         const IVec& arr,
                         const DVec& original_weight,
                         const OriginalEdgeIndex& original_edges,
                         Rcpp::NumericMatrix values,
                         Rcpp::NumericMatrix result)
    : graph(graph),
      dep(dep),
      arr(arr),
      original_weight(original_weight),
      original_edges(original_edges),
      values(values),
      result(result) {}

  int original_edge_index(int from, int to) const {
    int edge_index = -1;
    double best_weight = numeric_limits<double>::max();

    for (const auto& candidate : original_edges[from]) {
      int idx = candidate.second;
      if (candidate.first == to && original_weight[idx] < best_weight) {
        best_weight = original_weight[idx];
        edge_index = idx;
      }
    }

    return edge_index;
  }

  void operator()(std::size_t begin, std::size_t end) {
    DVec distances(graph->nbnode, numeric_limits<double>::max());
    DVec distances2(graph->nbnode, numeric_limits<double>::max());
    IVec visited1(graph->nbnode, 0);
    IVec visited2(graph->nbnode, 0);
    IVec visited;
    IVec parents(graph->nbnode, -1);
    IVec parents2(graph->nbnode, -1);

    for (std::size_t k = begin; k != end; ++k) {
      int start = dep[k];
      int target = arr[k];
      distances[start] = 0.0;
      distances2[target] = 0.0;
      PQ q;
      PQ qr;
      q.push(make_pair(start, 0.0));
      qr.push(make_pair(target, 0.0));
      int mid = -1;
      double mu = numeric_limits<double>::max();

      while (true) {
        bool q_done = q.empty() || q.top().second > mu;
        bool qr_done = qr.empty() || qr.top().second > mu;
        if (q_done && qr_done) break;
        if (q.empty() && qr.empty()) break;

        if (!q.empty()) {
          int v = q.top().first;
          double w = q.top().second;
          q.pop();

          visited.push_back(v);
          visited1[v] = 1;

          if (visited2[v] == 1 && distances[v] + distances2[v] < mu) {
            mid = v;
            mu = distances[v] + distances2[v];
          }

          if (w <= distances[v] && !Stall_par(v, distances, graph->nodeGr, graph->wGr, graph->indGr)) {
            for (int i = graph->indG[v]; i < graph->indG[v + 1]; ++i) {
              int v2 = graph->nodeG[i];
              double w2 = graph->wG[i];

              if (distances[v] + w2 < distances[v2]) {
                distances[v2] = distances[v] + w2;
                q.push(make_pair(v2, distances[v2]));
                visited1[v2] = 1;
                parents[v2] = v;
                visited.push_back(v2);
              }
            }
          }
        }

        if (!qr.empty()) {
          int v = qr.top().first;
          double w = qr.top().second;
          qr.pop();

          visited.push_back(v);
          visited2[v] = 1;

          if (visited1[v] == 1 && distances[v] + distances2[v] < mu) {
            mid = v;
            mu = distances[v] + distances2[v];
          }

          if (w <= distances2[v] && !Stall_par(v, distances2, graph->nodeG, graph->wG, graph->indG)) {
            for (int i = graph->indGr[v]; i < graph->indGr[v + 1]; ++i) {
              int v2 = graph->nodeGr[i];
              double w2 = graph->wGr[i];

              if (distances2[v] + w2 < distances2[v2]) {
                distances2[v2] = distances2[v] + w2;
                qr.push(make_pair(v2, distances2[v2]));
                visited2[v2] = 1;
                parents2[v2] = v;
                visited.push_back(v2);
              }
            }
          }
        }
      }

      if (mid != -1 && mu < numeric_limits<double>::max()) {
        IVec path;
        for (int p = parents2[mid]; p != -1; p = parents2[p]) {
          path.insert(path.begin(), p);
        }
        path.push_back(mid);
        for (int p = parents[mid]; p != -1; p = parents[p]) {
          path.push_back(p);
        }
        reverse(path.begin(), path.end());

        if (path.size() > 1) {
          graph->unpack(path);
        }

        result(k, 0) = mu;
        for (std::size_t j = 0; j < values.ncol(); ++j) {
          result(k, j + 1) = 0.0;
        }
        for (int i = 0; i < static_cast<int>(path.size()) - 1; ++i) {
          int edge_index = original_edge_index(path[i], path[i + 1]);
          if (edge_index == -1) continue;

          for (std::size_t j = 0; j < values.ncol(); ++j) {
            result(k, j + 1) += values(edge_index, j);
          }
        }
      }

      for (int node : visited) {
        visited1[node] = 0;
        visited2[node] = 0;
        distances[node] = numeric_limits<double>::max();
        distances2[node] = numeric_limits<double>::max();
        parents[node] = -1;
        parents2[node] = -1;
      }
      visited.clear();
    }
  }
};

// [[Rcpp::export]]
Rcpp::NumericMatrix cpppathvaluesC(std::vector<int> &orfrom,
                                   std::vector<int> &orto,
                                   std::vector<double> &orw,
                                   Rcpp::NumericMatrix values,
                                   std::vector<int> &gfrom,
                                   std::vector<int> &gto,
                                   std::vector<double> &gw,
                                   int nb,
                                   std::vector<int> &rank,
                                   std::vector<int> &shortf,
                                   std::vector<int> &shortt,
                                   std::vector<int> &shortc,
                                   bool phast,
                                   std::vector<int> dep,
                                   std::vector<int> arr) {
  CGraph network(gfrom, gto, gw, nb, rank, shortf, shortt, shortc, phast);
  network.construct_shortcuts();
  network.to_adj_list(false, phast);
  network.to_adj_list(true, phast);

  Rcpp::NumericMatrix result(dep.size(), values.ncol() + 1);
  std::fill(result.begin(), result.end(), Rcpp::NumericVector::get_na());

  OriginalEdgeIndex original_edges(nb);
  for (int i = 0; i < static_cast<int>(orfrom.size()); ++i) {
    original_edges[orfrom[i]].push_back(make_pair(orto[i], i));
  }

  CHPathValuesPairWorker worker(&network, dep, arr, orw,
                                original_edges, values, result);
  parallelFor(0, dep.size(), worker);

  return result;
}
