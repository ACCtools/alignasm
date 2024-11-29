//
// Created by pentagon03 on 8/12/24.
//

#ifndef ALIGNASM_GRAPH_OPERATIONS_HPP
#define ALIGNASM_GRAPH_OPERATIONS_HPP
#include <vector>

template<typename Dist_t>
using Graph = std::vector<std::vector<std::pair<int64_t, Dist_t>>>; // consider using vector<basic_string<Edge>> for performance boost

template<typename Dist_t>
void add_edge(Graph<Dist_t> &graph, int64_t from, int64_t to, Dist_t dist){
    assert(0<=from and from<graph.size() and 0<=to and to<graph.size());
    graph[from].push_back({to, dist});
}

template<typename Dist_t>
int64_t get_edge_count(const Graph<Dist_t>& graph){
    int64_t ans = 0;
    for(int64_t i = 0; i < graph.size(); i++)
        ans += graph[i].size();
    return ans;
}

#endif //ALIGNASM_GRAPH_OPERATIONS_HPP
