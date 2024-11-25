//
// Created by pentagon03 on 7/18/24.
//

#ifndef ALIGNASM_K_WEIGHTED_BFS_HPP
#define ALIGNASM_K_WEIGHTED_BFS_HPP

#include<vector>
#include<string>
#include<algorithm>
// Source: Dial's Algorithm
// link: https://codeforces.com/blog/entry/88408?#comment-790017
// backtest: https://www.acmicpc.net/source/share/27912f0d18284f28b3de1288e1d38d71

template<typename Graph>
void k_weighted_bfs(Graph& graph, int64_t src, int64_t lim, std::vector<int64_t> &dist, std::vector<int64_t> &pre) {
    ++lim; //[0...lim)
    std::vector<std::basic_string<int64_t>> qs(lim);
    dist.assign(graph.size(), -1);
    pre.assign(graph.size(), -1);

    dist[src] = 0; qs[0].push_back(src);
    for (int64_t d = 0, maxd = 0; d <= maxd; ++d) {
        for (auto& q = qs[d % lim]; !q.empty(); ) {
            int64_t cur = q.back(); q.pop_back();
            if (dist[cur] != d) continue;
            for (auto [nxt, cost] : graph[cur]) {
                assert(cost < lim);
                auto nd = d + cost;
                if (dist[nxt] != -1 && dist[nxt] <= nd) continue;
                dist[nxt] = nd; pre[nxt] = cur;
                qs[nd % lim].push_back(nxt);
                maxd = std::max(maxd, nd);
            }
        }
    }
}

#endif //ALIGNASM_K_WEIGHTED_BFS_HPP
