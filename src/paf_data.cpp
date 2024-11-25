#include "paf_data.hpp"

#include "graph_operations.hpp"
#include "k_shortest_walks.hpp"
#include "k_weighted_bfs.hpp"
#include "priority_queue_vector.hpp"

#include <iostream>
#include <charconv>
#include <limits>


void get_overlap_range(PafReadData &paf_read_data, std::string_view cs_str) {
    auto cs_iter = cs_str.begin() + CS_TAG_START;
    decltype(cs_iter) cs_chk_iter;
//    std::cout << cs_str << std::endl;

    int64_t ref_step;
    int64_t ref_index;
    int64_t qry_index = paf_read_data.qry_str;

    // str is the real begin
    ref_index = paf_read_data.ref_str;
    ref_step = paf_read_data.aln_fwd ? 1 : -1;

    int64_t converted_num = -1;
    paf_read_data.ref_overlap_range = std::vector<std::pair<int64_t, int64_t>> {};
    paf_read_data.qry_overlap_range = std::vector<std::pair<int64_t, int64_t>> {};
    while (cs_iter != cs_str.end()) {
        // Analysis CS tag
        // :7287-at:1234+t:242-tag:123*tg*ga

        if (*cs_iter == ':') {
            assert(cs_iter != cs_str.end() and std::next(cs_iter) != cs_str.end());
            cs_chk_iter = ++cs_iter;
            // Returns first non-digit pointer
            auto result_iter = std::from_chars(cs_iter, cs_str.end(), converted_num);
            cs_iter = result_iter.ptr;

            assert(cs_iter > cs_chk_iter);
            if (cs_iter > cs_chk_iter) {
                paf_read_data.ref_overlap_range.emplace_back(ref_index, ref_index + (converted_num - 1) * ref_step);
                ref_index += converted_num * ref_step;

                paf_read_data.qry_overlap_range.emplace_back(qry_index, qry_index + (converted_num - 1));
                qry_index += converted_num;
            }
        } else {
            cs_chk_iter = cs_iter++;
            while (cs_iter != cs_str.end() and isalpha(*cs_iter) > 0) {
                ++cs_iter;
            }

            auto variant_len = cs_iter - cs_chk_iter - 1;
            assert(variant_len > 0);

            if (*cs_chk_iter == '+') {
                qry_index += variant_len;
            } else if (*cs_chk_iter == '-') {
                ref_index += variant_len * ref_step;
            } else if (*cs_chk_iter == '*') {
                assert(variant_len == 2);
                ref_index += 1 * ref_step;
                qry_index += 1;
            } else {
                // Never happens
                assert(false);
            }
        }
    }

}



/// Get Best path in paf_ctg_data_original
void solve_ctg_read(std::vector<PafReadData> &paf_ctg_data_original, std::vector<PafOutputData> &paf_ctg_out, std::vector<PafOutputData> &paf_ctg_alt_out, std::vector<std::vector<PafOutputData>> &paf_ctg_max_out) {

    /// Test Sesson
    // do some tests here

    /// Start
    static int64_t test_case = 0;
    ++test_case;

    /// PAF_DATA SORT
    auto paf_ctg_data_sorted = paf_ctg_data_original;

    assert(not paf_ctg_data_sorted.empty()); // "There should at least be 1 data"
    if (paf_ctg_data_sorted.size() == 1) {
        paf_ctg_data_original[0].ctg_sorted_index = 0;
        paf_ctg_out.emplace_back(paf_ctg_data_sorted[0]); //emplace_back
        return;
    }

    std::sort(paf_ctg_data_sorted.begin(), paf_ctg_data_sorted.end()); // sort by (qry, ref)
    for(int64_t i = 0; i < (int64_t) paf_ctg_data_sorted.size(); i++){
        assert(0 <= paf_ctg_data_sorted[i].ctg_index and paf_ctg_data_sorted[i].ctg_index < paf_ctg_data_original.size());
        paf_ctg_data_original[paf_ctg_data_sorted[i].ctg_index].ctg_sorted_index = i;
    }
    // never do smth on paf_ctg_data_original

    /// paf_ctg_part: denotes continously overlapped ranges
    auto n = (int64_t) paf_ctg_data_sorted.size();
    std::vector<int32_t> paf_ctg_part{};
    std::vector<int64_t> part_idx(n);
    int64_t part_end = -1;
    for (int64_t idx = 0; idx < n; idx++) {
        const auto &read_data = paf_ctg_data_sorted[idx];
        if (part_end < read_data.qry_str) {  // doesn't overlap
            paf_ctg_part.push_back(idx);
        }
        part_idx[idx] = (int64_t) paf_ctg_part.size() - 1; // directly go to part index
        part_end = std::max(read_data.qry_end, part_end);
    }
    paf_ctg_part.push_back((int32_t) paf_ctg_data_sorted.size());

    /// Get Edited Loc after cutting
    // pair when failed to edit
    constexpr std::pair<int64_t, int64_t> FAIL_CUT{-1, -1};
    // (x -> y)
    // qry, ref for x in (x,y)
    std::vector<std::vector<std::pair<int64_t, int64_t>>> edited_loc_pre_end(n, std::vector<std::pair<int64_t, int64_t>>(n, FAIL_CUT));
    // qry, ref for y in (x,y)
    std::vector<std::vector<std::pair<int64_t, int64_t>>> edited_loc_str(n, std::vector<std::pair<int64_t, int64_t>>(n, FAIL_CUT));
    // (x qry index, y qry index)
    std::vector<std::vector<std::pair<int64_t, int64_t>>> edited_overlap_idx(n, std::vector<std::pair<int64_t, int64_t>>(n, FAIL_CUT));

    for (int64_t i = 0; i < n; i++) {
        assert(paf_ctg_data_sorted[i].qry_overlap_range.size() == paf_ctg_data_sorted[i].ref_overlap_range.size());
        for (const auto &[l, r]: paf_ctg_data_sorted[i].qry_overlap_range) {
            assert(l <= r);
        }
    }

    /// Vertex Index
    std::vector<std::vector<int64_t>> vtx_index(n, std::vector<int64_t>(n, -1));
    std::vector<std::pair<int64_t, int64_t>> index_vtx;

    // single one
    for (int64_t i = 0; i < n; i++) {
        vtx_index[i][i] = (int64_t) index_vtx.size();
        index_vtx.emplace_back(i, i);
        edited_loc_str[i][i] = {paf_ctg_data_sorted[i].qry_str, paf_ctg_data_sorted[i].ref_str};
        edited_overlap_idx[i][i] = {0, 0};
    }

    for (int64_t i = 0; i < n; i++) {
        const auto &pre = paf_ctg_data_sorted[i];
        int64_t pre_len = (int64_t) pre.qry_overlap_range.size();
        for (int64_t j = i + 1; j < n; j++) {
            const auto &cur = paf_ctg_data_sorted[j];
            if (pre.qry_end < cur.qry_str) break;
            assert(part_idx[i] == part_idx[j]);
            int64_t cur_len = (int64_t) cur.qry_overlap_range.size();
            if (qry_partial_overlap(pre, cur)) { // anom_pre doesn't contain cur
                bool determined = false;
                int64_t min_gap = -1;
                std::pair<int64_t, int64_t> min_gap_idx = {-1, -1};
                int64_t ref_step = cur.aln_fwd ? 1 : -1;
                int64_t ref_step_pre = pre.aln_fwd ? 1 : -1;
                for (int64_t p_i = 0, p_j = 0; p_i < pre_len and p_j < cur_len;) {
                    const auto &[l_i, r_i] = pre.qry_overlap_range[p_i];
                    const auto &[l_j, r_j] = cur.qry_overlap_range[p_j];
                    if (l_i == l_j) {
                        if (l_j == r_j) {
                            p_j++;
                            continue;
                        }
                        edited_loc_pre_end[i][j].first = l_i;
                        edited_loc_pre_end[i][j].second = pre.ref_overlap_range[p_i].first + (l_i - l_i) * ref_step_pre;
                        edited_loc_str[i][j].first = l_j + 1;
                        edited_loc_str[i][j].second = cur.ref_overlap_range[p_j].first + ((l_j + 1) - l_j) * ref_step;
                        edited_overlap_idx[i][j] = {p_i, p_j};
                        determined = true;
                        break;
                    }
                    if (l_i < l_j) {
                        if (l_j <= r_i + 1) {
                            edited_loc_pre_end[i][j].first = l_j - 1;
                            edited_loc_pre_end[i][j].second =
                                    pre.ref_overlap_range[p_i].first + ((l_j - 1) - l_i) * ref_step_pre;
                            edited_loc_str[i][j].first = l_j;
                            edited_loc_str[i][j].second = cur.ref_overlap_range[p_j].first + (l_j - l_j) * ref_step;
                            edited_overlap_idx[i][j] = {p_i, p_j};
                            determined = true;
                            break;
                        } else {
                            auto gap = l_j - (r_i + 1);
                            assert(gap >= 1);
                            if (min_gap == -1 or gap < min_gap) {
                                min_gap = gap;
                                min_gap_idx = {p_i, p_j};
                            }
                        }
                        p_i++;
                    } else {
                        if (l_i <= r_j - 1) {
                            edited_loc_pre_end[i][j].first = l_i;
                            edited_loc_pre_end[i][j].second = pre.ref_overlap_range[p_i].first + (l_i - l_i) * ref_step_pre;
                            edited_loc_str[i][j].first = l_i + 1;
                            edited_loc_str[i][j].second = cur.ref_overlap_range[p_j].first + (l_i + 1 - l_j) * ref_step;
                            edited_overlap_idx[i][j] = {p_i, p_j};
                            determined = true;
                            break;
                        }
                        p_j++;
                    }
                }
                if (determined or min_gap != -1) {
                    if (not determined) {
                        auto [p_i, p_j] = min_gap_idx;
                        const auto &[l_i, r_i] = pre.qry_overlap_range[p_i];
                        const auto &[l_j, r_j] = cur.qry_overlap_range[p_j];
                        edited_loc_pre_end[i][j].first = r_i;
                        edited_loc_pre_end[i][j].second = pre.ref_overlap_range[p_i].first + (r_i - l_i) * ref_step_pre;
                        edited_loc_str[i][j].first = l_j;
                        edited_loc_str[i][j].second = cur.ref_overlap_range[p_j].first + (l_j - l_j) * ref_step;
                        edited_overlap_idx[i][j] = {p_i, p_j};
                    }
                    vtx_index[i][j] = (int64_t) index_vtx.size();
                    index_vtx.emplace_back(i, j);
                } else {
                    assert(0 && "found a case that overlapping ith paf and jth paf can't be connected");
                }
            }
        }
    }

    /// Graph Construction
    int64_t vtx_n = (int64_t) index_vtx.size();
    auto vtx_to_index = [&](int64_t i, int64_t j) -> int64_t {
        assert(0 <= i and i < n and 0 <= j and j < n);
        return vtx_index[i][j];
    };
    auto index_to_vtx = [&](int64_t i) -> std::pair<int64_t, int64_t> {
        assert(0 <= i and i < index_vtx.size());
        return index_vtx[i];
    };

    // This stores edited qry range, ref range
    struct Internal_Vertex {
        int64_t pre_idx{-1};
        int64_t cur_idx{-1};
        bool is_one{true}; // is one idx
        int64_t qry_str{}, qry_end{};
        int64_t ref_str{}, ref_end{};
        bool default_vertex{false};
        Internal_Vertex() = default;
        void set_idx(int64_t idx){
            pre_idx = cur_idx = idx;
        }
        Internal_Vertex(int64_t i, int64_t j,
                        const std::vector<std::vector<std::pair<int64_t, int64_t>>> &edited_loc,
        const std::vector<PafReadData> &paf_data)
        : pre_idx(i), cur_idx(j), is_one(i == j),
        qry_str(edited_loc[i][j].first), ref_str(edited_loc[i][j].second),
        qry_end(paf_data[j].qry_end), ref_end(paf_data[j].ref_end), default_vertex(true){
            assert(0 <= i and i <= j and j < paf_data.size());
        }
    };
    auto linkable = [&](Internal_Vertex lft, Internal_Vertex rht) -> bool {
        assert(lft.default_vertex == rht.default_vertex);
        // we will just chk qrys, that's all
        if(not lft.default_vertex){
            return lft.qry_end < rht.qry_str;
        }

        assert(lft.default_vertex);
        // not valid vertices
        if (edited_loc_str[lft.pre_idx][lft.cur_idx] == FAIL_CUT or edited_loc_str[rht.pre_idx][rht.cur_idx] == FAIL_CUT)
            return false;
        if (not rht.is_one) {
            if (lft.cur_idx != rht.pre_idx) return false;
            return lft.qry_str < rht.qry_str;
        } else {
            if (part_idx[lft.cur_idx] + 1 == part_idx[rht.cur_idx]) return true;
            if (part_idx[lft.cur_idx] != part_idx[rht.cur_idx]) return false;
            return lft.qry_end < rht.qry_str; //no overlap
        }
    };

    // always call linkable min_distance
    // diff score => number of empty spaces between "lft" and "rht"
    // ex) diff([1,3] , [4,6]) = 4 - 3 - 1 = 0 */
    auto get_score = [&](Internal_Vertex lft, Internal_Vertex rht, bool calc_sum_) -> Paf_Distance {
        assert(linkable(lft, rht));
        auto ref_abs = [](auto x){
            if(x < 0){
                return -x * REF_NEGATIVE_PENALTY;
            }else{
                return x;
            }
        };
        Paf_Distance dist(calc_sum_);
        auto &[calc_sum, qry_score, ref_score, anom, qul_nonzero, qul_total] = dist;
        if (not rht.is_one) {
            assert(lft.cur_idx == rht.pre_idx);
            // let's change lft
            lft.qry_end = edited_loc_pre_end[rht.pre_idx][rht.cur_idx].first;
            lft.ref_end = edited_loc_pre_end[rht.pre_idx][rht.cur_idx].second;
        }
        // No Overlap Guranteed
        assert(lft.qry_end < rht.qry_str);
        int64_t qry_diff = rht.qry_str - lft.qry_end - 1;
        // we don't need SV_BASELINE
//        if (qry_diff > SV_BASELINE) {
//            anom += 1;
//            qry_diff = SV_BASELINE;
//        }
        int64_t ref_diff = 0;
        if (paf_ctg_data_sorted[lft.cur_idx].ref_chr == paf_ctg_data_sorted[rht.cur_idx].ref_chr
            and paf_ctg_data_sorted[lft.cur_idx].aln_fwd == paf_ctg_data_sorted[rht.cur_idx].aln_fwd) {
            // Case: same reference and same align direction
            if (paf_ctg_data_sorted[lft.cur_idx].aln_fwd) {
                // Case: forward direction
                ref_diff += ref_abs((rht.ref_str) - (lft.ref_end + 1));
            } else {
                // Case: Backward Direction
                ref_diff += ref_abs((rht.ref_str + 1) - (lft.ref_end));
            }
            if (ref_diff > SV_BASELINE) {
                anom += 1;
                ref_diff = SV_BASELINE;
            }
        } else if (paf_ctg_data_sorted[lft.cur_idx].ref_chr == paf_ctg_data_sorted[rht.cur_idx].ref_chr and
                   paf_ctg_data_sorted[lft.cur_idx].aln_fwd != paf_ctg_data_sorted[rht.cur_idx].aln_fwd) {
            // Case: different align direction
            anom += 1;
            ref_diff += SV_INV_PENALTY;
            if (paf_ctg_data_sorted[lft.cur_idx].aln_fwd) {
                ref_diff += ref_abs(rht.ref_end - (lft.ref_end + 1));
            } else {
                ref_diff += ref_abs(rht.ref_str - (lft.ref_str + 1));
            }
            if (ref_diff > SV_BASELINE) {
                anom += 1;
                ref_diff = SV_BASELINE;
            }
        } else {
            // Case: different reference OR
            assert(paf_ctg_data_sorted[lft.cur_idx].ref_chr != paf_ctg_data_sorted[rht.cur_idx].ref_chr);
            anom += 1;
            ref_diff = SV_TRANS_PENALTY;
        }
        assert(qry_diff >= 0 and ref_diff >= 0);
        qry_score = qry_diff * QRY_WEIGHT;
        ref_score = ref_diff * REF_WEIGHT;
        if (paf_ctg_data_sorted[rht.cur_idx].map_qul) qul_nonzero += 1;
        qul_total += 1;
        return dist;
    };

    /* Make Graph (DAG)
    Vertex: {paf_data index}
    Paf_Distance: {int64_t: Score(a,b), int64_t: ANOM(a, b)}
    Edge: // vertex x -> Edge{vertex, distance}
    Dijkstra Vertex: {distance, Vertex x}
    Add Virtual Node st, en */
    auto make_Graph = [&](Graph<Paf_Distance> &graph, int64_t n, int64_t src, int64_t dest) -> void {
        graph.clear(); graph.resize(n);
        /**
         * src_20240922_pentagon03 -> min_distance group
         * last group -> dest
         * ith group -> (i+1)th group
         */
        assert(paf_ctg_part.size() >= 2);
        // start -> min_distance
        {
            // [l, r)
            auto l = paf_ctg_part[0], r = paf_ctg_part[1];
            assert(l == 0);
            graph[src].reserve(r - l);
            for (auto i = l; i < r; i++) {
                Paf_Distance dist(true);
                auto &[calc_sum, qry_score, ref_score, anom, qul_nonzero, qul_total] = dist;
                qry_score += paf_ctg_data_sorted[i].qry_str * SV_FRONT_END_COEFFICIENT;
//                if (qry_score > SV_BASELINE) {
//                    qry_score = SV_BASELINE;
//                    anom += 1;
//                }
                if (paf_ctg_data_sorted[i].map_qul) qul_nonzero += 1;
                qul_total += 1;
                add_edge(graph, src, vtx_to_index(i, i), dist);
            }
        }
        // last -> end
        {
            // [l, r)
            assert(paf_ctg_part.size() >= 2);
            auto l = paf_ctg_part[(int64_t)paf_ctg_part.size() - 2], r = paf_ctg_part[(int64_t)paf_ctg_part.size() - 1];
            assert(r == paf_ctg_data_sorted.size());
            for (auto i = l; i < r; i++) {
                Paf_Distance dist(true);
                auto &[calc_sum, qry_score, ref_score, anom, qul_nonzero, qul_total] = dist;
                calc_sum = true;
                qry_score += (paf_ctg_data_sorted[i].qry_total_length - paf_ctg_data_sorted[i].qry_end - 1) * SV_FRONT_END_COEFFICIENT;
//                if (qry_score > SV_BASELINE) {
//                    qry_score = SV_BASELINE;
//                    anom += 1;
//                }
                add_edge(graph, vtx_to_index(i, i), dest, dist);
            }
        }
        // ith block Inside
        {
            for (int64_t block = 0; block + 1 < paf_ctg_part.size(); block++) {
                // Process current block
                auto l = paf_ctg_part[block], r = paf_ctg_part[block + 1];
                for (int64_t i = l; i < r; i++) {
                    for (int64_t j = i + 1; j < r; j++) {
                        if (paf_ctg_data_sorted[i].qry_contains(paf_ctg_data_sorted[j])) continue;
                        if (paf_ctg_data_sorted[i].qry_end < paf_ctg_data_sorted[j].qry_str) {
                            // no overlap, (i,i) -> (j,j)
                            Internal_Vertex IV_ii(i, i, edited_loc_str, paf_ctg_data_sorted);
                            Internal_Vertex IV_jj(j, j, edited_loc_str, paf_ctg_data_sorted);
                            if (linkable(IV_ii, IV_jj)) {
                                add_edge(graph, vtx_to_index(i, i), vtx_to_index(j, j), get_score(IV_ii, IV_jj, true));
                            }
                        } else {
                            // (i,i) -> (i, j)
                            Internal_Vertex IV_ii(i, i, edited_loc_str, paf_ctg_data_sorted);
                            Internal_Vertex IV_ij(i, j, edited_loc_str, paf_ctg_data_sorted);
                            if (linkable(IV_ii, IV_ij)) {
                                add_edge(graph, vtx_to_index(i, i), vtx_to_index(i, j), get_score(IV_ii, IV_ij, true));
                            }
                            // (i, j) -> (k, k),  (i, j) -> (j, k)
                            for (int64_t k = j + 1; k < r; k++) {
                                Internal_Vertex IV_kk(k, k, edited_loc_str, paf_ctg_data_sorted);
                                if (linkable(IV_ij, IV_kk)) {
                                    add_edge(graph, vtx_to_index(i, j), vtx_to_index(k, k), get_score(IV_ij, IV_kk, true));
                                }
                                Internal_Vertex IV_jk(j, k, edited_loc_str, paf_ctg_data_sorted);
                                if (linkable(IV_ij, IV_jk)) {
                                    add_edge(graph, vtx_to_index(i, j), vtx_to_index(j, k), get_score(IV_ij, IV_jk, true));
                                }
                            }
                        }
                    }
                }
            }
        }
        // ith block -> i+1 th block
        {
            for (int64_t block = 0; block + 2 < paf_ctg_part.size(); block++) {
                auto l = paf_ctg_part[block], r = paf_ctg_part[block + 1];
                auto l2 = paf_ctg_part[block + 1], r2 = paf_ctg_part[block + 2];
                for (int64_t i = l; i < r; i++) {
                    // (i, i) -> (k, k)
                    Internal_Vertex IV_ii(i, i, edited_loc_str, paf_ctg_data_sorted);
                    for (int64_t k = l2; k < r2; k++) {
                        Internal_Vertex IV_kk(k, k, edited_loc_str, paf_ctg_data_sorted);
                        if (linkable(IV_ii, IV_kk)) {
                            add_edge(graph, vtx_to_index(i, i), vtx_to_index(k, k), get_score(IV_ii, IV_kk,true));
                        }
                    }
                    // (i, j) -> (k, k)
                    for (int64_t j = i + 1; j < r; j++) {
                        if (paf_ctg_data_sorted[i].qry_contains(paf_ctg_data_sorted[j])) continue;
                        if (paf_ctg_data_sorted[i].qry_end < paf_ctg_data_sorted[j].qry_str) break;
                        Internal_Vertex IV_ij(i, j, edited_loc_str, paf_ctg_data_sorted);
                        for (int64_t k = l2; k < r2; k++) {
                            Internal_Vertex IV_kk(k, k, edited_loc_str, paf_ctg_data_sorted);
                            if (linkable(IV_ij, IV_kk)) {
                                add_edge(graph, vtx_to_index(i, j), vtx_to_index(k, k), get_score(IV_ij, IV_kk, true));
                            }
                        }
                    }
                }
            }
        }
    };

    /// Graph Construction
    int64_t src = vtx_n++;
    int64_t dest = vtx_n++;
    Graph<Paf_Distance> graph;
    make_Graph(graph, vtx_n, src, dest);

    /// ANOM GRAPH Construction
    Graph<int64_t> anom_graph(vtx_n);
    for (int64_t cur = 0; cur < vtx_n; cur++) {
        for (const auto &[nxt, dist]: graph[cur]) {
            add_edge<int64_t>(anom_graph, cur, nxt, dist.anom);
        }
    }
    constexpr int64_t MAX_ANOM = 2;
    std::vector<int64_t> anom_dis, anom_pre;
    k_weighted_bfs(anom_graph, src, MAX_ANOM, anom_dis, anom_pre);

    assert(anom_dis[dest] != -1);
    {
        // DO SMTH with anom_path
        int64_t cur = dest;
        std::vector<int64_t> anom_path;
        while (cur != -1) {
            anom_path.push_back(cur);
            cur = anom_pre[cur];
        }
        std::reverse(anom_path.begin(), anom_path.end());
    }

    /// Actual SubSequence of Pafs that has good (score, anom)
    kShortestWalksSolver sol(graph, Paf_Distance::max(), Paf_Distance(true), true, false);
    const int64_t MAX_PATH_COUNT = 10000; // Maximum Paths to watch for shorter anom score paths
    auto k_path_distances = sol.k_shortest_paths(src, dest, MAX_PATH_COUNT);

    assert(not k_path_distances.empty());
    auto min_distance = *k_path_distances.begin();

    /// Get Actual Sequence of Paf Subsequences that has good (score, anom)
    // there is path with (x,y) pairs Edges.
    // we deconstruct them with (x,x) or (x,y) and case work to make it a Pair with {x} vertices.
    using EdgePath = std::vector<std::tuple<int64_t, int64_t, Paf_Distance>>;
    using PafPath = std::vector<PafOutputData>;

    // Take note of Obsidian [[alignasm Project#7 31]]
    auto edge_path_to_paf_path = [&](const EdgePath &path) -> PafPath{
        assert(path.size() >= 2);
        PafPath paf_path{};
        for(const auto&[u, v, w]: path){
            if(u == src){
                auto [x, y] = index_to_vtx(v);
                assert(x>=0 and x < paf_ctg_data_sorted.size() and x == y);
                paf_path.emplace_back(paf_ctg_data_sorted[x]);
            }else if(v == dest){
                // do nothing
            }else{
                auto [x1, x2] = index_to_vtx(u);
                if(x1 == x2){
                    auto[y1, y2] = index_to_vtx(v);
                    if(y1 == y2){
                        int64_t x = x1, y = y1;
                        paf_path.emplace_back(paf_ctg_data_sorted[y]);
                    }else{
                        assert(x2 == y1);
                        int64_t x = y1, y = y2;
                        paf_path.emplace_back(paf_ctg_data_sorted[y]);
                        // paf_path shouldn't be invalidated
                        {
                            assert(paf_path.size() >= 2);
                            auto &px = paf_path[paf_path.size() - 2];
                            px.edited_qry_end = edited_loc_pre_end[x][y].first;
                            px.edited_ref_end = edited_loc_pre_end[x][y].second;
                            px.edited_qry_end_overlap_idx = edited_overlap_idx[x][y].first;
                            auto &py = paf_path[paf_path.size() - 1];
                            py.edited_qry_str = edited_loc_str[x][y].first;
                            py.edited_ref_str = edited_loc_str[x][y].second;
                            py.edited_qry_str_overlap_idx = edited_overlap_idx[x][y].second;
                        }
                    }
                }else{
                    auto [y1, y2] = index_to_vtx(v);
                    if(y1 == y2){
                        assert(x2 != y2);
                        int64_t z = y2;
                        paf_path.emplace_back(paf_ctg_data_sorted[z]);
                    }else{
                        int64_t x = x1, y = x2;
                        assert(y == y1);
                        int64_t z = y2;
                        paf_path.emplace_back(paf_ctg_data_sorted[z]);
                        assert(paf_path.size() >= 2);
                        // paf_path shouldn't be invalidated
                        {
                            auto &py = paf_path[paf_path.size() - 2];
                            py.edited_qry_end = edited_loc_pre_end[y][z].first;
                            py.edited_ref_end = edited_loc_pre_end[y][z].second;
                            py.edited_qry_end_overlap_idx = edited_overlap_idx[y][z].first;
                            auto &pz = paf_path[paf_path.size() - 1];
                            pz.edited_qry_str = edited_loc_str[y][z].first;
                            pz.edited_ref_str = edited_loc_str[y][z].second;
                            pz.edited_qry_str_overlap_idx = edited_overlap_idx[y][z].second;
                        }
                    }
                }
            }
        }
        return paf_path;
    };
    auto InternalVertex_to_PafOutputData = [&](const Internal_Vertex& IV){
        PafOutputData node;
        node.ctg_index = paf_ctg_data_sorted[IV.cur_idx].ctg_index;
        node.edited_qry_str = IV.qry_str;
        node.edited_qry_end = IV.qry_end;
        node.edited_ref_str = IV.ref_str;
        node.edited_ref_end = IV.ref_end;
        return node;
    };
    // We find the gaps and fill it with larger piece
    auto upgrade_paf_path = [&](const PafPath &paf_path) -> PafPath{
        int64_t len = (int64_t) paf_path.size();
        assert(len >= 1);
        auto upgraded_paf_path = PafPath{};
        int64_t qry_min = std::numeric_limits<int64_t>::max(), qry_max = std::numeric_limits<int64_t>::min();
        for(const auto &read_data: paf_ctg_data_sorted) {
            qry_min = std::min(qry_min, read_data.qry_str);
            qry_max = std::max(qry_max, read_data.qry_end);
        }
        PQVec<std::pair<int64_t, int64_t>, std::vector<std::pair<int64_t, int64_t> >, std::greater<> > pq;
        int64_t sorted_idx_iter = 0;
        auto add_first_node = [&](){
            auto cur = paf_path[0];
            auto l = (int64_t) qry_min, r = cur.edited_qry_str - 1;
            if(l > r){
                upgraded_paf_path.push_back(cur);
                return;
            }
            if(l == r){
                upgraded_paf_path.push_back(cur);
                return;
            }
            while(not pq.empty() && pq.top().first < r)
                pq.pop();
            while(sorted_idx_iter < (int64_t) paf_ctg_data_sorted.size() && paf_ctg_data_sorted[sorted_idx_iter].qry_str <= l){
                if(paf_ctg_data_sorted[sorted_idx_iter].qry_end >= r)
                    pq.emplace(paf_ctg_data_sorted[sorted_idx_iter].qry_end, sorted_idx_iter);
                sorted_idx_iter++;
            }
            if(pq.empty()){
                upgraded_paf_path.push_back(cur);
                return;
            }
            const auto &rdata = paf_ctg_data_original[cur.ctg_index];
            const auto& candidates = pq.getVector();
            Paf_Distance min_dist = Paf_Distance::max();
            Internal_Vertex ans_IV, ans_IVR;
            bool ans_flag = false;
            for(auto [end_, sorted_idx]: candidates){
                assert(end_ >= r);
                // Apply Internal Vertex
                Internal_Vertex IV, IVR;

                IV.set_idx(sorted_idx);

                IVR.set_idx(rdata.ctg_sorted_index);
                IVR.qry_end = cur.edited_qry_end;
                IVR.ref_end = cur.edited_ref_end;
                const auto &data = paf_ctg_data_sorted[sorted_idx];
                const auto &qry_range = data.qry_overlap_range;
                bool ok_l = false, ok_r = false;
                int64_t range_idx = 0;
                // l
                {
                    assert(l == data.qry_str);
                    IV.qry_str = data.qry_str;
                    IV.ref_str = data.ref_str;
                    ok_l = true;
                }
                // r
                {
                    auto it = std::lower_bound(qry_range.begin(), qry_range.end(), std::pair<int64_t, int64_t>(r+1, -1));
                    if(it != qry_range.begin() && prev(it)->second >= r){
                        it = prev(it);
                        assert(it->first <= r and r <= it->second);
                        int64_t nxt_range_idx = (int64_t)(it - qry_range.begin());
                        assert(r >= l);
                        assert(nxt_range_idx >= range_idx);
                        range_idx = nxt_range_idx;
                        IVR.qry_str = cur.edited_qry_str;
                        IVR.qry_end = cur.edited_qry_end;
                        IVR.ref_str = cur.edited_ref_str;
                        IVR.ref_end = cur.edited_ref_end;
                        IV.qry_end = r;
                        int64_t ref_step = data.aln_fwd ? 1 : -1;
                        assert(0 <= range_idx and range_idx < data.ref_overlap_range.size());
                        IV.ref_end =  data.ref_overlap_range[range_idx].first + (r - it->first) * ref_step;
                        ok_r = true;
                    }else{
                        bool determined = false;
                        const auto &r_qry_range = rdata.qry_overlap_range;
                        int64_t ref_step_pre = data.aln_fwd ? 1 : -1;
                        int64_t ref_step = rdata.aln_fwd ? 1 : -1;
                        std::pair<int64_t, int64_t> edited_loc_pre_end, edited_loc_str;
                        for(int64_t p = 0; p < (int64_t) r_qry_range.size() and range_idx < (int64_t)qry_range.size(); ){
                            if(r_qry_range[p].first > cur.edited_qry_end){
                                break;
                            }
                            auto [l_i, r_i] = qry_range[range_idx];
                            auto [l_j, r_j] = r_qry_range[p];
                            if(r_j > cur.edited_qry_end){
                                r_j = cur.edited_qry_end;
                            }
                            if(l_i == l_j){
                                if(l_j == r_j){
                                    range_idx++;
                                    continue;
                                }
                                edited_loc_pre_end.first = l_i;
                                edited_loc_pre_end.second = data.ref_overlap_range[range_idx].first + (l_i - l_i) * ref_step_pre;
                                edited_loc_str.first = l_j + 1;
                                edited_loc_str.second = rdata.ref_overlap_range[p].first + ((l_j + 1) - l_j) * ref_step;
                                determined = true;
                                break;
                            }
                            if (l_i < l_j) {
                                if (l_j <= r_i + 1) {
                                    edited_loc_pre_end.first = l_j - 1;
                                    edited_loc_pre_end.second =
                                            data.ref_overlap_range[range_idx].first + ((l_j - 1) - l_i) * ref_step_pre;
                                    edited_loc_str.first = l_j;
                                    edited_loc_str.second = rdata.ref_overlap_range[p].first + (l_j - l_j) * ref_step;
                                    determined = true;
                                    break;
                                }
                                range_idx++;
                            } else {
                                if (l_i <= r_j - 1) {
                                    edited_loc_pre_end.first = l_i;
                                    edited_loc_pre_end.second =
                                            data.ref_overlap_range[range_idx].first + (l_i - l_i) * ref_step_pre;
                                    edited_loc_str.first = l_i + 1;
                                    edited_loc_str.second =
                                            rdata.ref_overlap_range[p].first + (l_i + 1 - l_j) * ref_step;
                                    determined = true;
                                    break;
                                }
                                p++;
                            }
                        }
                        if(determined){
                            IV.qry_end = edited_loc_pre_end.first;
                            IV.ref_end = edited_loc_pre_end.second;
                            IVR.qry_str = edited_loc_str.first;
                            IVR.ref_str = edited_loc_str.second;
                            ok_r = true;
                        }
                    }
                }
                if(not ok_r) continue;
                assert(ok_l and ok_r);
                assert(not IV.default_vertex and not IVR.default_vertex);
                auto score = get_score(IV, IVR, true);
                if(score < min_dist){
                    min_dist = score;
                    ans_IV = IV;
                    ans_IVR = IVR;
                    ans_flag = true;
                }
            }
            if(not ans_flag){
                upgraded_paf_path.push_back(cur);
                return;
            }else{
                assert(upgraded_paf_path.empty());
                PafOutputData pre_node, cur_node, nxt_node;
                cur_node = InternalVertex_to_PafOutputData(ans_IV);
                nxt_node = InternalVertex_to_PafOutputData(ans_IVR);
                upgraded_paf_path.push_back(cur_node);
                upgraded_paf_path.push_back(nxt_node);
            }
        };
        auto add_main_nodes = [&](){
            // process (i -> i+1) nodes
            for(int64_t i = 1; i < len; i++){
                // pre contains suffix
                auto pre = upgraded_paf_path.back();
                // cur contains prefix
                auto cur = paf_path[i];
                auto l = pre.edited_qry_end + 1, r = cur.edited_qry_str - 1;
                if(l > r) {
                    upgraded_paf_path.push_back(cur);
                    continue;
                }
                if(l == r){
                    upgraded_paf_path.push_back(cur);
                    continue;
                }
                while(not pq.empty() && pq.top().first < r)
                    pq.pop();
                while(sorted_idx_iter < (int64_t) paf_ctg_data_sorted.size() && paf_ctg_data_sorted[sorted_idx_iter].qry_str <= l){
                    if(paf_ctg_data_sorted[sorted_idx_iter].qry_end >= r)
                        pq.emplace(paf_ctg_data_sorted[sorted_idx_iter].qry_end, sorted_idx_iter);
                    sorted_idx_iter++;
                }
                if(pq.empty()){
                    upgraded_paf_path.push_back(cur);
                    continue;
                }
                const auto &ldata = paf_ctg_data_original[pre.ctg_index];
                const auto &rdata = paf_ctg_data_original[cur.ctg_index];
                const auto& candidates = pq.getVector();
                Paf_Distance min_dist = Paf_Distance::max();
                Internal_Vertex ans_IVL, ans_IV, ans_IVR;
                bool ans_flag = false;
                for(auto [end_, sorted_idx]: candidates){
                    assert(end_ >= r);
                    // Apply Internal Vertex
                    Internal_Vertex IVL, IV, IVR;

                    IVL.set_idx(ldata.ctg_sorted_index);
                    IVL.qry_str = pre.edited_qry_str;
                    IVL.ref_str = pre.edited_ref_str;

                    IV.set_idx(sorted_idx);

                    IVR.set_idx(rdata.ctg_sorted_index);
                    IVR.qry_end = cur.edited_qry_end;
                    IVR.ref_end = cur.edited_ref_end;
                    const auto &data = paf_ctg_data_sorted[sorted_idx];
                    const auto &qry_range = data.qry_overlap_range;
                    bool ok_l = false, ok_r = false;
                    int64_t range_idx = 0;

                    // l
                    {
                        auto it = std::lower_bound(qry_range.begin(), qry_range.end(), std::pair<int64_t, int64_t>(l+1, -1));
                        if(it != qry_range.begin() && prev(it)->second >= l){
                            it = prev(it);
                            assert(it->first <= l and l <= it->second);
                            range_idx = (int64_t)(it - qry_range.begin());
                            IVL.qry_str = pre.edited_qry_str;
                            IVL.qry_end = pre.edited_qry_end;
                            IVL.ref_str = pre.edited_ref_str;
                            IVL.ref_end = pre.edited_ref_end;
                            IV.qry_str = l;
                            int64_t ref_step = data.aln_fwd ? 1 : -1;
                            assert(0 <= range_idx and range_idx < data.ref_overlap_range.size());
                            IV.ref_str =  data.ref_overlap_range[range_idx].first + (l - it->first) * ref_step;
                            ok_l = true;
                        }else{
                            bool determined = false;
                            const auto &l_qry_range = ldata.qry_overlap_range;
                            int64_t ref_step_pre = ldata.aln_fwd ? 1 : -1;
                            int64_t ref_step = data.aln_fwd ? 1 : -1;
                            std::pair<int64_t, int64_t> edited_loc_pre_end, edited_loc_str;
                            for(int64_t p = 0; p < (int64_t) l_qry_range.size() and range_idx < (int64_t)qry_range.size(); ){
                                if(l_qry_range[p].second < pre.edited_qry_str){
                                    p++;
                                    continue;
                                }
                                auto [l_i, r_i] = l_qry_range[p];
                                auto [l_j, r_j] = qry_range[range_idx];
                                if(l_i < pre.edited_qry_str){
                                    l_i = pre.edited_qry_str;
                                }
                                if(l_i == l_j){
                                    if(l_j == r_j){
                                        range_idx++;
                                        continue;
                                    }
                                    edited_loc_pre_end.first = l_i;
                                    edited_loc_pre_end.second = ldata.ref_overlap_range[p].first + (l_i - l_qry_range[p].first) * ref_step_pre;
                                    edited_loc_str.first = l_j + 1;
                                    edited_loc_str.second = data.ref_overlap_range[range_idx].first + ((l_j + 1) - l_j) * ref_step;
                                    determined = true;
                                    break;
                                }
                                if (l_i < l_j) {
                                    if (l_j <= r_i + 1) {
                                        edited_loc_pre_end.first = l_j - 1;
                                        edited_loc_pre_end.second =
                                                ldata.ref_overlap_range[p].first + ((l_j - 1) - l_qry_range[p].first) * ref_step_pre;
                                        edited_loc_str.first = l_j;
                                        edited_loc_str.second = data.ref_overlap_range[range_idx].first + (l_j - l_j) * ref_step;
                                        determined = true;
                                        break;
                                    }
                                    p++;
                                } else {
                                    if (l_i <= r_j - 1) {
                                        edited_loc_pre_end.first = l_i;
                                        edited_loc_pre_end.second =
                                                ldata.ref_overlap_range[p].first + (l_i - l_qry_range[p].first) * ref_step_pre;
                                        edited_loc_str.first = l_i + 1;
                                        edited_loc_str.second =
                                                data.ref_overlap_range[range_idx].first + (l_i + 1 - l_j) * ref_step;
                                        determined = true;
                                        break;
                                    }
                                    range_idx++;
                                }
                            }
                            if(determined){
                                IVL.qry_end = edited_loc_pre_end.first;
                                IVL.ref_end = edited_loc_pre_end.second;
                                IV.qry_str = edited_loc_str.first;
                                IV.ref_str = edited_loc_str.second;
                                ok_l = true;
                            }
                        }
                    }
                    if(not ok_l) continue;

                    // r
                    {
                        auto it = std::lower_bound(qry_range.begin(), qry_range.end(), std::pair<int64_t, int64_t>(r+1, -1));
                        if(it != qry_range.begin() && prev(it)->second >= r){
                            it = prev(it);
                            assert(it->first <= r and r <= it->second);
                            int64_t nxt_range_idx = (int64_t)(it - qry_range.begin());
                            assert(r >= l);
                            assert(nxt_range_idx >= range_idx);
                            range_idx = nxt_range_idx;
                            IVR.qry_str = cur.edited_qry_str;
                            IVR.qry_end = cur.edited_qry_end;
                            IVR.ref_str = cur.edited_ref_str;
                            IVR.ref_end = cur.edited_ref_end;
                            IV.qry_end = r;
                            int64_t ref_step = data.aln_fwd ? 1 : -1;
                            assert(0 <= range_idx and range_idx < data.ref_overlap_range.size());
                            IV.ref_end =  data.ref_overlap_range[range_idx].first + (r - it->first) * ref_step;
                            ok_r = true;
                        }else{
                            bool determined = false;
                            const auto &r_qry_range = rdata.qry_overlap_range;
                            int64_t ref_step_pre = data.aln_fwd ? 1 : -1;
                            int64_t ref_step = rdata.aln_fwd ? 1 : -1;
                            std::pair<int64_t, int64_t> edited_loc_pre_end, edited_loc_str;
                            for(int64_t p = 0; p < (int64_t) r_qry_range.size() and range_idx < (int64_t)qry_range.size(); ){
                                if(r_qry_range[p].first > cur.edited_qry_end){
                                    break;
                                }
                                auto [l_i, r_i] = qry_range[range_idx];
                                auto [l_j, r_j] = r_qry_range[p];
                                if(r_j > cur.edited_qry_end){
                                    r_j = cur.edited_qry_end;
                                }
                                if(l_i == l_j){
                                    if(l_j == r_j){
                                        range_idx++;
                                        continue;
                                    }
                                    edited_loc_pre_end.first = l_i;
                                    edited_loc_pre_end.second = data.ref_overlap_range[range_idx].first + (l_i - l_i) * ref_step_pre;
                                    edited_loc_str.first = l_j + 1;
                                    edited_loc_str.second = rdata.ref_overlap_range[p].first + ((l_j + 1) - l_j) * ref_step;
                                    determined = true;
                                    break;
                                }
                                if (l_i < l_j) {
                                    if (l_j <= r_i + 1) {
                                        edited_loc_pre_end.first = l_j - 1;
                                        edited_loc_pre_end.second =
                                                data.ref_overlap_range[range_idx].first + ((l_j - 1) - l_i) * ref_step_pre;
                                        edited_loc_str.first = l_j;
                                        edited_loc_str.second = rdata.ref_overlap_range[p].first + (l_j - l_j) * ref_step;
                                        determined = true;
                                        break;
                                    }
                                    range_idx++;
                                } else {
                                    if (l_i <= r_j - 1) {
                                        edited_loc_pre_end.first = l_i;
                                        edited_loc_pre_end.second =
                                                data.ref_overlap_range[range_idx].first + (l_i - l_i) * ref_step_pre;
                                        edited_loc_str.first = l_i + 1;
                                        edited_loc_str.second =
                                                rdata.ref_overlap_range[p].first + (l_i + 1 - l_j) * ref_step;
                                        determined = true;
                                        break;
                                    }
                                    p++;
                                }
                            }
                            if(determined){
                                IV.qry_end = edited_loc_pre_end.first;
                                IV.ref_end = edited_loc_pre_end.second;
                                IVR.qry_str = edited_loc_str.first;
                                IVR.ref_str = edited_loc_str.second;
                                ok_r = true;
                            }
                        }
                    }
                    if(not ok_r) continue;
                    assert(ok_l and ok_r);
                    assert(not IVL.default_vertex and not IV.default_vertex and not IVR.default_vertex);
                    auto score = get_score(IVL, IV, true) + get_score(IV, IVR, true);
                    if(score < min_dist){
                        min_dist = score;
                        ans_IVL = IVL;
                        ans_IV = IV;
                        ans_IVR = IVR;
                        ans_flag = true;
                    }
                }
                if(not ans_flag){
                    upgraded_paf_path.push_back(cur);
                    continue;
                }else{
                    assert(not upgraded_paf_path.empty());
                    PafOutputData pre_node, cur_node, nxt_node;
                    upgraded_paf_path.pop_back();
                    pre_node = InternalVertex_to_PafOutputData(ans_IVL);
                    cur_node = InternalVertex_to_PafOutputData(ans_IV);
                    nxt_node = InternalVertex_to_PafOutputData(ans_IVR);
                    upgraded_paf_path.push_back(pre_node);
                    upgraded_paf_path.push_back(cur_node);
                    upgraded_paf_path.push_back(nxt_node);
                }
            }
        };
        // add last node
        auto add_last_node = [&](){
            auto pre = upgraded_paf_path.back();
            auto l = pre.edited_qry_end + 1, r = qry_max;
            if(l > r){
                return;
            }
            if(l == r){
                return;
            }
            while(not pq.empty() && pq.top().first < r)
                pq.pop();
            while(sorted_idx_iter < (int64_t) paf_ctg_data_sorted.size() && paf_ctg_data_sorted[sorted_idx_iter].qry_str <= l){
                if(paf_ctg_data_sorted[sorted_idx_iter].qry_end >= r)
                    pq.emplace(paf_ctg_data_sorted[sorted_idx_iter].qry_end, sorted_idx_iter);
                sorted_idx_iter++;
            }
            if(pq.empty()){
                return;
            }
            const auto &ldata = paf_ctg_data_original[pre.ctg_index];
            const auto& candidates = pq.getVector();
            Paf_Distance min_dist = Paf_Distance::max();
            Internal_Vertex ans_IVL, ans_IV;
            bool ans_flag = false;
            for(auto [end_, sorted_idx]: candidates){
                assert(end_ >= r);
                // Apply Internal Vertex
                Internal_Vertex IVL, IV;

                IVL.set_idx(ldata.ctg_sorted_index);
                IVL.qry_str = pre.edited_qry_str;
                IVL.ref_str = pre.edited_ref_str;

                IV.set_idx(sorted_idx);

                const auto &data = paf_ctg_data_sorted[sorted_idx];
                const auto &qry_range = data.qry_overlap_range;
                bool ok_l = false, ok_r = false;
                int64_t range_idx = 0;
                // l
                {
                    auto it = std::lower_bound(qry_range.begin(), qry_range.end(), std::pair<int64_t, int64_t>(l+1, -1));
                    if(it != qry_range.begin() && prev(it)->second >= l){
                        it = prev(it);
                        assert(it->first <= l and l <= it->second);
                        range_idx = (int64_t)(it - qry_range.begin());
                        IVL.qry_str = pre.edited_qry_str;
                        IVL.qry_end = pre.edited_qry_end;
                        IVL.ref_str = pre.edited_ref_str;
                        IVL.ref_end = pre.edited_ref_end;
                        IV.qry_str = l;
                        int64_t ref_step = data.aln_fwd ? 1 : -1;
                        assert(0 <= range_idx and range_idx < data.ref_overlap_range.size());
                        IV.ref_str =  data.ref_overlap_range[range_idx].first + (l - it->first) * ref_step;
                        ok_l = true;
                    }else{
                        bool determined = false;
                        const auto &l_qry_range = ldata.qry_overlap_range;
                        int64_t ref_step_pre = ldata.aln_fwd ? 1 : -1;
                        int64_t ref_step = data.aln_fwd ? 1 : -1;
                        std::pair<int64_t, int64_t> edited_loc_pre_end, edited_loc_str;
                        for(int64_t p = 0; p < (int64_t) l_qry_range.size() and range_idx < (int64_t)qry_range.size(); ){
                            if(l_qry_range[p].second < pre.edited_qry_str){
                                p++;
                                continue;
                            }
                            auto [l_i, r_i] = l_qry_range[p];
                            auto [l_j, r_j] = qry_range[range_idx];
                            if(l_i < pre.edited_qry_str){
                                l_i = pre.edited_qry_str;
                            }
                            if(l_i == l_j){
                                if(l_j == r_j){
                                    range_idx++;
                                    continue;
                                }
                                edited_loc_pre_end.first = l_i;
                                edited_loc_pre_end.second = ldata.ref_overlap_range[p].first + (l_i - l_qry_range[p].first) * ref_step_pre;
                                edited_loc_str.first = l_j + 1;
                                edited_loc_str.second = data.ref_overlap_range[range_idx].first + ((l_j + 1) - l_j) * ref_step;
                                determined = true;
                                break;
                            }
                            if (l_i < l_j) {
                                if (l_j <= r_i + 1) {
                                    edited_loc_pre_end.first = l_j - 1;
                                    edited_loc_pre_end.second =
                                            ldata.ref_overlap_range[p].first + ((l_j - 1) - l_qry_range[p].first) * ref_step_pre;
                                    edited_loc_str.first = l_j;
                                    edited_loc_str.second = data.ref_overlap_range[range_idx].first + (l_j - l_j) * ref_step;
                                    determined = true;
                                    break;
                                }
                                p++;
                            } else {
                                if (l_i <= r_j - 1) {
                                    edited_loc_pre_end.first = l_i;
                                    edited_loc_pre_end.second =
                                            ldata.ref_overlap_range[p].first + (l_i - l_qry_range[p].first) * ref_step_pre;
                                    edited_loc_str.first = l_i + 1;
                                    edited_loc_str.second =
                                            data.ref_overlap_range[range_idx].first + (l_i + 1 - l_j) * ref_step;
                                    determined = true;
                                    break;
                                }
                                range_idx++;
                            }
                        }
                        if(determined){
                            IVL.qry_end = edited_loc_pre_end.first;
                            IVL.ref_end = edited_loc_pre_end.second;
                            IV.qry_str = edited_loc_str.first;
                            IV.ref_str = edited_loc_str.second;
                            ok_l = true;
                        }
                    }
                }
                if(not ok_l) continue;
                // r
                {
                    assert(r == data.qry_end);
                    IV.qry_end = data.qry_end;
                    IV.ref_end = data.ref_end;
                    ok_r = true;
                }
                assert(ok_l and ok_r);
                assert(not IVL.default_vertex and not IV.default_vertex);
                auto score = get_score(IVL, IV, true);
                if(score < min_dist){
                    min_dist = score;
                    ans_IVL = IVL;
                    ans_IV = IV;
                    ans_flag = true;
                }
            }
            if(not ans_flag){
                return;
            }else{
                assert(not upgraded_paf_path.empty());
                PafOutputData pre_node, cur_node, nxt_node;
                upgraded_paf_path.pop_back();
                pre_node = InternalVertex_to_PafOutputData(ans_IVL);
                cur_node = InternalVertex_to_PafOutputData(ans_IV);
                upgraded_paf_path.push_back(pre_node);
                upgraded_paf_path.push_back(cur_node);
            }
        };
        add_first_node();
        add_main_nodes();
        add_last_node();
        return upgraded_paf_path;
    };

    /// Find EdgePath 1
    auto path1 = sol.kth_shortest_full_path(src, dest, 0, false);
    auto paf_path1 = edge_path_to_paf_path(path1);
    auto upgraded_paf_path1 = upgrade_paf_path(paf_path1);

    paf_ctg_out = upgraded_paf_path1;

    /// Find Max Edge Paths
    {
        int64_t idx = 1;
        auto is_equal_paf_distance = [](const Paf_Distance &lft, const Paf_Distance &rht)->bool{
            return lft.score_sum() == rht.score_sum() and lft.anom == rht.anom;
        };
        for(;idx<k_path_distances.size() && is_equal_paf_distance(min_distance, k_path_distances[idx]);idx++){
            auto path_max = sol.kth_shortest_full_path(src, dest, idx, false);
            auto paf_path_max = edge_path_to_paf_path(path_max);
            auto upgraded_paf_path_max = upgrade_paf_path(paf_path_max);
            paf_ctg_max_out.push_back(upgraded_paf_path_max);
        }
    }

    /// Find EdgePath 2 (Alt)
    if ((int64_t) k_path_distances.size() >= 2 and min_distance.anom != anom_dis[dest]) {
        Paf_Distance ans{true, -1, -1};
        int64_t ans_up{}, ans_down{};
        int64_t ans_idx = -1;
        for (int64_t i = 1; i < k_path_distances.size(); i++) {
            const auto &d = k_path_distances[i];
            if (d.anom >= min_distance.anom) continue;
            auto up = d.score_sum() - min_distance.score_sum();
            assert(up > 0);
            auto down = min_distance.anom - d.anom;
            assert(down > 0);
            if (ans_idx == -1 or up * ans_down < down * ans_up) {
                ans = d;
                ans_up = up;
                ans_down = down;
                ans_idx = i;
            }
        }

        if(ans_idx != -1) {
//            std::cout << min_distance.score_sum() << ',' << ans.score_sum() << ',' << paf_ctg_data_original[0].paf_index << '\n';
            auto path2 = sol.kth_shortest_full_path(src, dest, ans_idx, false);
            auto paf_path2 = edge_path_to_paf_path(path2);
            auto upgraded_paf_path2 = upgrade_paf_path(paf_path2);
            paf_ctg_alt_out = upgraded_paf_path2;
//            for (auto& paf_line : paf_path2){
//                auto& paf_data_line = paf_ctg_data_original[paf_line.ctg_index];
//
//                auto out = std::vector<std::string>({std::to_string(paf_ctg_data_original[0].paf_index + 1), std::to_string(paf_data_line.qry_total_length), std::to_string(paf_line.edited_qry_str), std::to_string(paf_line.edited_qry_end + 1),
//                                                     paf_data_line.aln_fwd ? "+" : "-", std::to_string(paf_data_line.ref_chr),
//                                                     std::to_string(paf_data_line.ref_total_length), std::to_string(paf_line.edited_ref_str), std::to_string(paf_line.edited_ref_end + 1)});
//                for (auto& t: out) {
//                    std::cout << t << '\t';
//                }
//                std::cout << '\n';
//            }
        }
    }
}