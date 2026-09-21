#ifndef ALIGNASM_QUALITY_SCORING_HPP
#define ALIGNASM_QUALITY_SCORING_HPP

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

enum class ScoringMode { QUALITY, LEGACY };

struct ScoringConfig {
    ScoringMode mode = ScoringMode::QUALITY;
    int64_t sv_cost = 2000;
    int64_t mapq_loss_per_kb = 10;
    void validate() const;
};

constexpr int64_t QUALITY_SCORE_SCALE = 60000;
int64_t checked_score_integer(__int128 value);
int64_t effective_mapq(uint8_t mapq);

// Query intervals are half-open. Reference endpoints follow query orientation
// and are inclusive, matching the coordinates stored in PafReadData.
struct QualityInterval {
    int64_t start, end;
    uint8_t mapq;
};

class MapqEnvelope {
public:
    MapqEnvelope(int64_t query_length, const std::vector<QualityInterval>& intervals);
    int64_t integral(int64_t end) const;
    int64_t loss(int64_t start, int64_t end, uint8_t mapq = 0) const;
    int64_t query_length() const { return query_length_; }
private:
    int64_t query_length_;
    std::vector<int64_t> positions_, levels_, prefix_;
};

struct QualitySegment {
    int64_t query_start, query_end;
    int64_t ref_first, ref_last;
    int32_t chromosome;
    bool forward;
    uint8_t mapq;
};

// All ranking coordinates are additive. In particular, there is no ratio or
// rounding in the ordering used by the shortest-path/sidetrack algorithms.
struct QualityDistance {
    int64_t total_scaled = 0;
    int64_t mapq_loss_scaled = 0;
    int64_t unmapped_bp = 0;
    int64_t anom = 0;
    int64_t pieces = 0;
    int64_t query_cost = 0;
    int64_t sv_cost = 0;
    int64_t mapq_deficit = 0;
    bool infinite = false;

    static QualityDistance max();
    auto key() const {
        return std::tie(total_scaled, mapq_loss_scaled, unmapped_bp, anom, pieces);
    }
    bool operator==(const QualityDistance& other) const;
    bool operator!=(const QualityDistance& other) const { return !(*this == other); }
    bool operator<(const QualityDistance& other) const;
    bool operator>(const QualityDistance& other) const { return other < *this; }
    bool operator<=(const QualityDistance& other) const { return !(other < *this); }
    bool operator>=(const QualityDistance& other) const { return !(*this < other); }
    QualityDistance operator+(const QualityDistance& other) const;
    QualityDistance operator-(const QualityDistance& other) const;
    bool same_components(const QualityDistance& other) const;
};

class QualityScorer {
public:
    QualityScorer(const ScoringConfig& config, int64_t query_length,
                  const std::vector<QualityInterval>& intervals);
    QualityDistance start(int64_t first_query_start) const;
    QualityDistance transition(const QualitySegment& left, const QualitySegment& right) const;
    QualityDistance finish(const QualitySegment& last) const;
    // Independent whole-path accounting, used to verify recovered graph paths.
    QualityDistance evaluate(const std::vector<QualitySegment>& path) const;
private:
    ScoringConfig config_;
    MapqEnvelope envelope_;
    void validate_segment(const QualitySegment& segment) const;
    std::pair<int64_t, int64_t> junction(const QualitySegment& left,
                                       const QualitySegment& right) const;
    QualityDistance distance(int64_t query_cost, int64_t sv_cost, int64_t deficit,
                             int64_t unmapped, int64_t joins, int64_t pieces) const;
};

struct QualityPathReport {
    std::string kind;
    int64_t path_number;
    QualityDistance score;
};

struct QualityReport {
    int64_t paths_examined = 0;
    bool limit_reached = false;
    std::vector<QualityPathReport> paths;
};

#endif
