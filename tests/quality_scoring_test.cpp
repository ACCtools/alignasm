#include "quality_scoring.hpp"
#include "graph_operations.hpp"
#include "k_shortest_walks.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>

namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

template<typename Exception, typename Function>
void throws(Function function, const char* message) {
    try { function(); }
    catch (const Exception&) { return; }
    throw std::runtime_error(message);
}

QualitySegment piece(int64_t start, int64_t end, int64_t reference, uint8_t mapq = 60,
                     int32_t chromosome = 0, bool forward = true) {
    const auto last = reference + end - start - 1;
    return {start, end, forward ? reference : last, forward ? last : reference,
            chromosome, forward, mapq};
}

void test_envelope() {
    MapqEnvelope envelope(20, {{2, 12, 7}, {5, 10, 60}, {8, 15, 255}, {12, 16, 80}});
    require(envelope.integral(0) == 0, "empty integral");
    require(envelope.integral(5) == 21, "prefix integral");
    require(envelope.integral(10) == 321, "overlap maximum");
    require(envelope.integral(12) == 335, "ending high-quality interval");
    require(envelope.integral(20) == 575, "capped MAPQ and missing suffix");
    require(envelope.loss(5, 10) == 300, "unmapped interval loses evidence");
    require(envelope.loss(2, 5, 7) == 0, "retained evidence");
    require(envelope.loss(12, 16, 80) == 0, "MAPQ above 60 is capped");
    require(effective_mapq(255) == 0 && effective_mapq(254) == 60, "MAPQ normalization");
    MapqEnvelope unknown(100, {{0, 100, 255}, {0, 100, 0}});
    require(unknown.loss(0, 100, 255) == 0, "unknown quality is neutral");
    MapqEnvelope long_query(1000000000000LL, {});
    require(long_query.integral(1000000000000LL) == 0, "sparse profile for long query");
    throws<std::logic_error>([&] { envelope.loss(0, 1, 60); }, "invalid selected quality must fail");
    throws<std::invalid_argument>([&] { envelope.loss(2, 1); }, "reversed query interval must fail");
    throws<std::invalid_argument>([] { MapqEnvelope invalid(10, {{0, 11, 60}}); }, "out-of-range candidate must fail");
}

void test_real_coordinates_and_omission() {
    const std::vector<QualityInterval> candidates{{0, 23843, 7}, {23847, 48647, 60},
                                                 {0, 23847, 0}, {23847, 48647, 0}};
    QualityScorer scorer({}, 48647, candidates);
    const std::vector<QualitySegment> primary{
        {0, 23843, 148634036, 148657866, 0, true, 7},
        {23847, 48647, 83127712, 83102912, 0, false, 60}};
    const std::vector<QualitySegment> old_selection{
        {0, 23847, 145213613, 145188421, 0, false, 0},
        {23847, 48647, 78031620, 78006377, 1, false, 0}};
    const auto good = scorer.evaluate(primary);
    const auto bad = scorer.evaluate(old_selection);
    require(good.total_scaled == 2004 * QUALITY_SCORE_SCALE, "primary cost must be 2004");
    require(good.mapq_deficit == 0 && good.unmapped_bp == 4 && good.anom == 1, "primary components");
    require(bad.total_scaled == 136549010 && bad.mapq_deficit == 1654901, "old-selection analytical score");
    require(good < bad, "four bp must not defeat long high-quality anchors");
    for (const auto& path : {primary, old_selection}) {
        const auto graph_score = scorer.start(path.front().query_start) +
            scorer.transition(path[0], path[1]) + scorer.finish(path.back());
        require(graph_score.same_components(scorer.evaluate(path)), "edge/path accounting agreement");
    }
    const auto omitted = scorer.evaluate({primary.front()});
    require(omitted.mapq_deficit == 24800 * 60 && omitted.mapq_loss_scaled == 14880000,
            "omitting an anchor must retain its MAPQ penalty");
    require(omitted.unmapped_bp == 24804 && omitted.query_cost == 49608, "terminal omission cost");
}

void test_clipping_splitting_and_zero_weight() {
    QualityScorer scorer({}, 100, {{0, 100, 60}, {0, 100, 0}});
    auto clipped = piece(30, 70, 100, 60);
    const auto score = scorer.evaluate({clipped});
    require(score.mapq_deficit == 3600 && score.query_cost == 120 && score.unmapped_bp == 60,
            "clipping must lose prefix/suffix quality exactly once");
    require(score.same_components(scorer.start(30) + scorer.finish(clipped)), "clipped graph verification");
    const auto full = scorer.evaluate({piece(0, 100, 100)});
    const auto split = scorer.evaluate({piece(0, 40, 100), piece(40, 100, 140)});
    require(full.total_scaled == split.total_scaled && full.mapq_deficit == split.mapq_deficit,
            "splitting cannot create quality rewards");
    require(full < split, "fewer pieces break a genuine tie");

    QualityScorer disabled({ScoringMode::QUALITY, 2000, 0}, 100, {{0, 100, 60}, {0, 100, 0}});
    require(disabled.evaluate({piece(0, 100, 100, 60)}) == disabled.evaluate({piece(0, 100, 100, 0)}),
            "K=0 disables the quality tie-break as well as the main penalty");

    // One bp at MAPQ 7 has a fractional cost; it must survive until the final sum.
    QualityScorer fractional({}, 2, {{0, 2, 7}, {0, 2, 0}});
    const auto fractional_split = fractional.evaluate({piece(0, 1, 100, 0), piece(1, 2, 101, 0)});
    require(fractional_split.mapq_loss_scaled == 140 && fractional_split.total_scaled == 140,
            "never round a per-piece quality penalty");
}

void test_junctions_and_symmetry() {
    const std::vector<QualityInterval> candidates{{0, 10, 60}, {20, 30, 60}, {0, 30, 0}};
    QualityScorer scorer({}, 30, candidates);
    const auto left = piece(0, 10, 100);
    require(scorer.transition(left, piece(20, 30, 120)).anom == 0, "matching query/ref gaps are collinear");
    require(scorer.transition(left, piece(20, 30, 10)).sv_cost == 100, "negative gap is not doubled");
    require(scorer.transition(left, piece(20, 30, 10000000)).sv_cost == 2000, "large gap saturates");
    const auto inversion = scorer.transition(left, piece(20, 30, 10000000, 60, 0, false));
    require(inversion.sv_cost == 2000 && inversion.anom == 1, "inversion is one capped junction");
    const auto translocation = scorer.transition(left, piece(20, 30, 5, 60, 1));
    require(translocation.sv_cost == 2000 && translocation.anom == 1, "translocation is neutral to chromosome number");
    for (const auto right : {piece(20, 30, 10), piece(20, 30, 120),
                            piece(20, 30, 10000000, 60, 0, false), piece(20, 30, 5, 60, 1)}) {
        const std::vector<QualitySegment> path{left, right};
        auto shifted = path;
        for (auto& node : shifted) { node.ref_first += 9000; node.ref_last += 9000; }
        require(scorer.evaluate(path).same_components(scorer.evaluate(shifted)), "reference translation symmetry");
        auto reversed = path;
        std::reverse(reversed.begin(), reversed.end());
        for (auto& node : reversed) {
            const auto start = node.query_start;
            node.query_start = 30 - node.query_end;
            node.query_end = 30 - start;
            std::swap(node.ref_first, node.ref_last);
            node.forward = !node.forward;
        }
        require(scorer.evaluate(path).same_components(scorer.evaluate(reversed)), "reverse-complement symmetry");
    }
}

void test_checked_arithmetic() {
    const auto maximum = std::numeric_limits<int64_t>::max();
    throws<std::overflow_error>([&] { MapqEnvelope envelope(maximum, {{0, maximum, 60}}); }, "envelope overflow");
    throws<std::overflow_error>([&] {
        QualityScorer scorer({ScoringMode::QUALITY, maximum, 10}, 20, {{0, 20, 60}});
        scorer.transition(piece(0, 10, 0), piece(10, 20, 0, 60, 1));
    }, "SV score overflow");
    throws<std::overflow_error>([&] {
        QualityScorer scorer({ScoringMode::QUALITY, 2000, maximum}, 20, {{0, 20, 60}});
        scorer.finish(piece(0, 20, 0, 0));
    }, "MAPQ score overflow");
    QualityDistance large;
    large.total_scaled = maximum;
    QualityDistance one;
    one.total_scaled = 1;
    throws<std::overflow_error>([&] { auto unused = large + one; }, "distance addition overflow");
    require(one < QualityDistance::max() && !(QualityDistance::max() < one), "infinity ordering");
    require(((large - one) + one).same_components(large), "signed sidetrack arithmetic");
}

void test_exhaustive_dag_oracle() {
    std::mt19937 random(20260921);
    for (int trial = 0; trial < 150; ++trial) {
        const int count = 2 + random() % 7;
        Graph<QualityDistance> graph(count);
        for (int from = 0; from < count; ++from) {
            for (int to = from + 1; to < count; ++to) {
                if (to != from + 1 && random() % 3) continue;
                QualityDistance weight;
                weight.query_cost = random() % 3;
                weight.sv_cost = random() % 3;
                weight.mapq_deficit = random() % 5;
                weight.mapq_loss_scaled = (trial % 2 ? 10 : 0) * weight.mapq_deficit;
                weight.unmapped_bp = weight.query_cost;
                weight.anom = random() % 2;
                weight.pieces = 1;
                weight.total_scaled = QUALITY_SCORE_SCALE * (weight.query_cost + weight.sv_cost) + weight.mapq_loss_scaled;
                add_edge(graph, from, to, weight);
            }
        }
        std::vector<QualityDistance> oracle;
        std::function<void(int, QualityDistance)> visit = [&](int node, QualityDistance distance) {
            if (node == count - 1) { oracle.push_back(distance); return; }
            for (const auto& [to, weight] : graph[node]) visit(static_cast<int>(to), distance + weight);
        };
        visit(0, {});
        std::sort(oracle.begin(), oracle.end());
        kShortestWalksSolver solver(graph, QualityDistance::max(), QualityDistance{}, true, false, true);
        const auto actual = solver.k_shortest_walks(0, count - 1, 10000);
        require(actual == oracle, "k-shortest distances must match exhaustive enumeration");
        for (size_t rank = 0; rank < actual.size(); ++rank) {
            const auto path = solver.kth_shortest_walk_recover(0, count - 1, rank);
            QualityDistance verified;
            int64_t previous = 0;
            for (const auto& [from, to, weight] : path) {
                require(from == previous, "recovery must be connected");
                const auto edge = std::find_if(graph[from].begin(), graph[from].end(), [&](const auto& item) { return item.first == to; });
                require(edge != graph[from].end() && weight.same_components(edge->second), "recovered edge must have original components");
                verified = verified + edge->second;
                previous = to;
            }
            require(previous == count - 1 && verified.same_components(actual[rank]), "recovered path must sum to its distance");
        }
        kShortestWalksSolver fast(graph, QualityDistance::max(), QualityDistance{}, true, false, true);
        require(fast.k_shortest_walks(0, count - 1, 1).front() == oracle.front(), "fast optimum is exact");
        require(fast.kth_shortest_walk_recover(0, count - 1, 0) == solver.kth_shortest_walk_recover(0, count - 1, 0),
                "--write-all cannot change the selected tied path");
    }
}
}

int main() {
    try {
        test_envelope();
        test_real_coordinates_and_omission();
        test_clipping_splitting_and_zero_weight();
        test_junctions_and_symmetry();
        test_checked_arithmetic();
        test_exhaustive_dag_oracle();
        std::cout << "quality scoring and exhaustive DAG checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
