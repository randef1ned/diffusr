/** diffusr: network diffusion algorithms in R
 *
 * Copyright (C) 2016 Simon Dirmeier
 * @author Simon Dirmeier
 * @email simon.dirmeier@bsse.ethz.ch
 *
 * This file is part of diffusr.
 *
 * diffusr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * diffusr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with diffusr. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../inst/include/diffusr.h"

struct distance_comparator
{
    bool operator()(pair<int, double>& lhs, pair<int, double>& rhs)
    {
        return lhs.second < rhs.second;
    }
};

bool equals(const double val, const double cmp, const double delta)
{
    return val <= cmp + delta && val >= cmp - delta;
}

priority_queue<pair<int, double>,
               vector<pair<int, double>>,
               distance_comparator>
  init_queue(const NumericMatrix& W, const int source)
{
    priority_queue<pair<int, double>,
                   vector<pair<int, double>>,
                   distance_comparator> queue;

    for (int i = 0; i < W.cols(); ++i)
    {
        if (i != source && W(source, i) > 0)
        {
            queue.push(make_pair(i, W(source, i)));
        }
    }

    return queue;
}

vector<pair<int, double>> current_neighbors(
  priority_queue<pair<int, double>,
                 vector<pair<int, double>>,
                 distance_comparator>& queue,
  vector<uint8_t>& visited)
{

    // get the nearest neighbor
    pair<int, double> cn = queue.top();
    // remove from queue
    queue.pop();

    vector<pair<int, double>> curr_nei;
    // add the node itself to the list, since we iterate over the list later
    if (!visited[cn.first])
    {
        curr_nei.push_back(cn);
    }

    // add neighbors that are as close as the first neighbor `cn`
    while (queue.size() && equals(queue.top().second, cn.second, .001))
    {
        pair<int, double> nn = queue.top();
        if (!visited[nn.first])
            curr_nei.push_back(nn);
        queue.pop();
    }

    return curr_nei;
}

void add_neighbor_to_queue(
  priority_queue<pair<int, double>,
                 vector<pair<int, double>>,
                 distance_comparator>& queue,
  const NumericMatrix& W,
  const pair<int, double>& cn)
{
    for (int i = 0; i < W.cols(); ++i)
    {
        if (i != cn.first && W(cn.first, i) > 0)
        {
            queue.push(make_pair(i, W(cn.first, i)));
        }
    }
}

set<int> nearest_neighbor_dijkstra_(
                                const int            source,
                                const int            max_depth,
                                const NumericMatrix& W)
{
    // use a priority queue to quickly extract nearest neighbors
    priority_queue<pair<int, double>,
                   vector<pair<int, double>>,
                   distance_comparator>
      queue = init_queue(W, source);
    // boolean vector if nodes have already been visited in the BFS
    vector<uint8_t> visited(W.rows(), false);

    // traverse graph until a certain depth is reached
    int r = 1;
    set<int> nei;
    do
    {
        checkUserInterrupt();
        // list of the nearest neighbors
        vector<pair<int, double>> curr_nei =
          current_neighbors(queue, visited);
        // iterate over current nearest neighbors and add them to the results
        // list
        for (const pair<int, double>& cn : curr_nei)
        {
            if (visited[cn.first] && cn.first == source)
            {
                continue;
            }
            else
            {
                visited[cn.first] = true;
            }
            // add as neighbor with index + 1, cause neighbors are indexe
            // starting from 1 in R
            nei.insert(cn.first + 1);
            // add current node to priority queue
            add_neighbor_to_queue(queue, W, cn);
        }
    } while (r++ < max_depth && queue.size());
    return nei;
}

//' Find the closest neighbors of a group of nodes in a graph.
//'
//' @noRd
//' @param node_idxs  the staring distribution
//' @param W  adjacency matrix
//' @param k  the depth of the nearest neighbor search
//' @return  returns a list of nearest neighbors for every node idxs given in
//'  <emph>node_idxs</emph>
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
List neighbors_(const vector<int> &node_idxs,
                const NumericMatrix& W,
                const int            k)
{
    // number of idxs given
    size_t len = node_idxs.size();
    // neighbors for every node
    vector<set<int>> neighbors(len);

    // parallelize node search
    // Fucking R bug: you can never have R objects within the parallel sections
    // #pragma omp parallel for
    for (uint32_t i = 0; i < len; ++i)
    {
        // substract one, cause R was one-based
        int node_idx = node_idxs[i] - 1;
        // neighbors of current node
        // run disjkstra until k neighbors are found
        neighbors[i] = nearest_neighbor_dijkstra_(node_idx, k, W);
    }

    return wrap(neighbors);
}
