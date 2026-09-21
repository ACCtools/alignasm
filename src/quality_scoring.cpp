#include "quality_scoring.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

void ScoringConfig::validate() const {
    if (sv_cost <= 0) throw std::invalid_argument("--sv-cost must be positive");
    if (mapq_loss_per_kb < 0)
        throw std::invalid_argument("--mapq-loss-per-kb must be nonnegative");
}

int64_t checked_score_integer(__int128 value) {
    if (value < std::numeric_limits<int64_t>::min() ||
        value > std::numeric_limits<int64_t>::max())
        throw std::overflow_error("quality score exceeds the signed 64-bit integer range");
    return static_cast<int64_t>(value);
}

int64_t effective_mapq(uint8_t mapq) {
    return mapq == 255 ? 0 : std::min<int64_t>(mapq, 60);
}

MapqEnvelope::MapqEnvelope(int64_t query_length,
                           const std::vector<QualityInterval>& intervals)
    : query_length_(query_length) {
    if (query_length <= 0) throw std::invalid_argument("query length must be positive");
    using Event = std::tuple<int64_t, int64_t, int64_t>;
    std::vector<Event> events{{0, 0, 0}, {query_length, 0, 0}};
    for (const auto& interval : intervals) {
        if (interval.start < 0 || interval.start >= interval.end || interval.end > query_length)
            throw std::invalid_argument("alignment query interval is outside the query");
        const auto quality = effective_mapq(interval.mapq);
        if (quality) {
            events.emplace_back(interval.start, quality, 1);
            events.emplace_back(interval.end, quality, -1);
        }
    }
    std::sort(events.begin(), events.end());
    std::array<int64_t, 61> counts{};
    int64_t previous = 0, level = 0, area = 0;
    for (size_t index = 0; index < events.size();) {
        const auto position = std::get<0>(events[index]);
        area = checked_score_integer(static_cast<__int128>(area) +
                                      static_cast<__int128>(position - previous) * level);
        while (index < events.size() && std::get<0>(events[index]) == position) {
            const auto [unused, quality, delta] = events[index++];
            counts[quality] += delta;
        }
        level = 60;
        while (level && counts[level] == 0) --level;
        positions_.push_back(position);
        levels_.push_back(level);
        prefix_.push_back(area);
        previous = position;
    }
}

int64_t MapqEnvelope::integral(int64_t end) const {
    if (end < 0 || end > query_length_)
        throw std::invalid_argument("MAPQ integral endpoint is outside the query");
    const auto index = static_cast<size_t>(
        std::upper_bound(positions_.begin(), positions_.end(), end) - positions_.begin() - 1);
    return checked_score_integer(static_cast<__int128>(prefix_[index]) +
        static_cast<__int128>(end - positions_[index]) * levels_[index]);
}

int64_t MapqEnvelope::loss(int64_t start, int64_t end, uint8_t mapq) const {
    if (start < 0 || end < start || end > query_length_)
        throw std::invalid_argument("invalid MAPQ loss interval");
    const auto value = checked_score_integer(static_cast<__int128>(integral(end)) -
        integral(start) - static_cast<__int128>(end - start) * effective_mapq(mapq));
    if (value < 0) throw std::logic_error("selected MAPQ exceeds its candidate envelope");
    return value;
}

QualityDistance QualityDistance::max() {
    QualityDistance distance;
    distance.infinite = true;
    return distance;
}

bool QualityDistance::operator==(const QualityDistance& other) const {
    if (infinite || other.infinite) return infinite == other.infinite;
    return key() == other.key();
}

bool QualityDistance::operator<(const QualityDistance& other) const {
    if (infinite) return false;
    if (other.infinite) return true;
    return key() < other.key();
}

namespace {
constexpr std::array<int64_t QualityDistance::*, 8> components{
    &QualityDistance::total_scaled, &QualityDistance::mapq_loss_scaled,
    &QualityDistance::unmapped_bp, &QualityDistance::anom, &QualityDistance::pieces,
    &QualityDistance::query_cost, &QualityDistance::sv_cost, &QualityDistance::mapq_deficit
};

QualityDistance combine(const QualityDistance& left, const QualityDistance& right, int sign) {
    if (left.infinite || right.infinite)
        throw std::logic_error("arithmetic on an infinite quality distance");
    QualityDistance result;
    for (const auto member : components)
        result.*member = checked_score_integer(static_cast<__int128>(left.*member) +
                                               static_cast<__int128>(sign) * (right.*member));
    return result;
}
}

QualityDistance QualityDistance::operator+(const QualityDistance& other) const {
    return combine(*this, other, 1);
}

QualityDistance QualityDistance::operator-(const QualityDistance& other) const {
    return combine(*this, other, -1);
}

bool QualityDistance::same_components(const QualityDistance& other) const {
    if (infinite != other.infinite) return false;
    for (const auto member : components)
        if (this->*member != other.*member) return false;
    return true;
}

QualityScorer::QualityScorer(const ScoringConfig& config, int64_t query_length,
                           const std::vector<QualityInterval>& intervals)
    : config_(config), envelope_(query_length, intervals) {
    config_.validate();
}

void QualityScorer::validate_segment(const QualitySegment& segment) const {
    if (segment.query_start < 0 || segment.query_start >= segment.query_end ||
        segment.query_end > envelope_.query_length() || segment.ref_first < 0 ||
        segment.ref_last < 0 ||
        (segment.forward ? segment.ref_first > segment.ref_last : segment.ref_first < segment.ref_last))
        throw std::invalid_argument("invalid quality-scoring alignment segment");
}

std::pair<int64_t, int64_t> QualityScorer::junction(
    const QualitySegment& left, const QualitySegment& right) const {
    if (left.query_end > right.query_start)
        throw std::logic_error("recovered query alignments overlap");
    if (left.chromosome != right.chromosome || left.forward != right.forward)
        return {config_.sv_cost, 1};
    const __int128 reference_gap = left.forward
        ? static_cast<__int128>(right.ref_first) - left.ref_last - 1
        : static_cast<__int128>(left.ref_last) - right.ref_first - 1;
    const auto absolute_gap = reference_gap < 0 ? -reference_gap : reference_gap;
    const auto cost = static_cast<int64_t>(std::min<__int128>(absolute_gap, config_.sv_cost));
    return {cost, reference_gap != right.query_start - left.query_end ? 1 : 0};
}

QualityDistance QualityScorer::distance(int64_t query_cost, int64_t sv_cost, int64_t deficit,
                                       int64_t unmapped, int64_t joins, int64_t pieces) const {
    QualityDistance score;
    score.query_cost = query_cost;
    score.sv_cost = sv_cost;
    score.mapq_deficit = deficit;
    score.unmapped_bp = unmapped;
    score.anom = joins;
    score.pieces = pieces;
    score.mapq_loss_scaled = checked_score_integer(
        static_cast<__int128>(config_.mapq_loss_per_kb) * deficit);
    score.total_scaled = checked_score_integer(
        QUALITY_SCORE_SCALE * (static_cast<__int128>(query_cost) + sv_cost) + score.mapq_loss_scaled);
    return score;
}

QualityDistance QualityScorer::start(int64_t first_query_start) const {
    const auto deficit = envelope_.loss(0, first_query_start);
    return distance(checked_score_integer(static_cast<__int128>(first_query_start) * 2),
                    0, deficit, first_query_start, 0, 0);
}

QualityDistance QualityScorer::transition(const QualitySegment& left,
                                          const QualitySegment& right) const {
    validate_segment(left);
    validate_segment(right);
    const auto [sv_cost, joins] = junction(left, right);
    const auto gap = right.query_start - left.query_end;
    const auto deficit = checked_score_integer(
        static_cast<__int128>(envelope_.loss(left.query_start, left.query_end, left.mapq)) +
        envelope_.loss(left.query_end, right.query_start));
    return distance(gap, sv_cost, deficit, gap, joins, 1);
}

QualityDistance QualityScorer::finish(const QualitySegment& last) const {
    validate_segment(last);
    const auto suffix = envelope_.query_length() - last.query_end;
    const auto deficit = checked_score_integer(
        static_cast<__int128>(envelope_.loss(last.query_start, last.query_end, last.mapq)) +
        envelope_.loss(last.query_end, envelope_.query_length()));
    return distance(checked_score_integer(static_cast<__int128>(suffix) * 2),
                    0, deficit, suffix, 0, 1);
}

QualityDistance QualityScorer::evaluate(const std::vector<QualitySegment>& path) const {
    if (path.empty()) throw std::logic_error("cannot score an empty alignment path");
    __int128 covered = 0, retained_quality = 0, reference_cost = 0, joins = 0;
    for (size_t index = 0; index < path.size(); ++index) {
        const auto& segment = path[index];
        validate_segment(segment);
        const auto length = segment.query_end - segment.query_start;
        covered += length;
        retained_quality += static_cast<__int128>(length) * effective_mapq(segment.mapq);
        if (index) {
            const auto [cost, count] = junction(path[index - 1], segment);
            reference_cost += cost;
            joins += count;
        }
    }
    const auto unmapped = checked_score_integer(envelope_.query_length() - covered);
    const auto deficit = checked_score_integer(envelope_.integral(envelope_.query_length()) - retained_quality);
    if (unmapped < 0 || deficit < 0) throw std::logic_error("invalid recovered quality coverage");
    // Unmapped bases count once internally and twice at the two ends.
    const auto query_cost = checked_score_integer(static_cast<__int128>(unmapped) +
        path.front().query_start + envelope_.query_length() - path.back().query_end);
    return distance(query_cost, checked_score_integer(reference_cost), deficit, unmapped,
                    checked_score_integer(joins), checked_score_integer(path.size()));
}
