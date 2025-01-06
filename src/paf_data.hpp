#ifndef ALIGNASM_PAF_DATA_HPP
#define ALIGNASM_PAF_DATA_HPP

#include <cassert>
#include <cinttypes>
#include <vector>
#include <string>
#include <string_view>
#include <iostream>

constexpr bool NON_SKIP_LINKABLE = true;
enum class UpgradeMode{
    ALT_PATH,
    SINGLE_PIECE,
    NO_UPDATE,
};
constexpr UpgradeMode UPGRADE_MODE = UpgradeMode::ALT_PATH;

constexpr int64_t QRY_WEIGHT = 1;
constexpr int64_t REF_WEIGHT = 1;
constexpr int64_t REF_NEGATIVE_PENALTY = 2;

constexpr int32_t CS_TAG_START = 5;
constexpr int64_t SV_BASELINE = 1e6;
constexpr int64_t SV_TRANS_PENALTY = 2000;
constexpr int64_t SV_INV_PENALTY = 500;
constexpr int64_t SV_FRONT_END_COEFFICIENT = 2;

// constexpr int64_t QRY_END_INF =

enum PafIndex {
    PAF_QRY_CHR,
    PAF_QRY_TOT,
    PAF_QRY_STR,
    PAF_QRY_END,
    PAF_ALN_FWD,
    PAF_REF_CHR,
    PAF_REF_TOT,
    PAF_REF_STR,
    PAF_REF_END,
    PAF_MAT_NUM,
    PAF_ALN_LEN,
    PAF_MAT_QUL
};

struct PafReadData {
    int32_t paf_index;
    int32_t ctg_index;
    int32_t ctg_sorted_index;
    std::string cs_string;
    int32_t mat_num;
    int32_t aln_len;
    std::vector<std::pair<int64_t, int64_t>> ref_overlap_range; // [l, r]
    std::vector<std::pair<int64_t, int64_t>> qry_overlap_range; // [l, r]
    int64_t ref_total_length; // the maximum range of ref.
    int64_t qry_total_length; // the maximum range of qry.
    int64_t qry_str, qry_end; // [qry_str, qry_end]
    int64_t ref_str, ref_end; // [ref_str, ref_end]
    int32_t ref_chr; // reference id
    uint8_t map_qul; // map quality
    bool aln_fwd;
    bool operator <(const PafReadData& rht) const{
        if(qry_str != rht.qry_str)
            return qry_str < rht.qry_str;
        return qry_end < rht.qry_end;
    }
    [[nodiscard]] bool qry_contains(const PafReadData &rht) const{
        // lft contains rht
        return qry_str <= rht.qry_str and rht.qry_end <= qry_end;
    }
    friend bool qry_partial_overlap(const PafReadData &lft, const PafReadData &rht){
        if(lft.qry_str < rht.qry_str){
            return rht.qry_str <= lft.qry_end && lft.qry_end < rht.qry_end;
        }else if(rht.qry_str < lft.qry_str){
            return lft.qry_str <= rht.qry_end && rht.qry_end < lft.qry_end;
        }else{
            return false;
        }
    }
};


struct PafOutputData {
    int32_t ctg_index;
    int64_t edited_qry_str, edited_qry_end;
    int64_t edited_ref_str, edited_ref_end;
    bool is_alt_path;
    PafOutputData()
            : ctg_index(-1), edited_qry_str(0), edited_qry_end(0),
              edited_ref_str(0), edited_ref_end(0), is_alt_path(false) {}
    void set_alt_path(bool res = true){
        is_alt_path = res;
    }
    explicit PafOutputData(const PafReadData& readData):
            ctg_index(readData.ctg_index),
            edited_qry_str(readData.qry_str), edited_qry_end(readData.qry_end),
            edited_ref_str(readData.ref_str), edited_ref_end(readData.ref_end), is_alt_path(false){}
};


struct PafEditData {
    std::string edit_cs_string;
    int32_t mat_num;
    int32_t aln_len;
    bool is_cut;
};

enum class PafDistanceCompareMode{
    CALC_SUM_MODE,
    QRY_SCORE_MODE,
};

/// Distance between Paf nodes
struct PafDistance{
    static thread_local PafDistanceCompareMode cmp_mode;
    static void set_mode(PafDistanceCompareMode mode){
        cmp_mode = mode;
    }
    // somehow, shouldn't initialize these below to use basic_string
    bool calc_sum_chk;
    int64_t qry_score, ref_score;
    int64_t anom;
    int64_t qul_nonzero, qul_total;
    PafDistance()=default;
    explicit PafDistance(bool calc_sum_chk_ , int64_t qry_score_ = 0, int64_t ref_score_ = 0,
                         int64_t anom_ = 0, int64_t qul_nonzero_ = 0, int64_t qul_total_ = 0)
            : calc_sum_chk(calc_sum_chk_), qry_score(qry_score_), ref_score(ref_score_),
              anom(anom_), qul_nonzero(qul_nonzero_), qul_total(qul_total_){}
    [[nodiscard]] static PafDistance max(){
        return PafDistance(false, -1, -1, -1, -1);
    }
    [[nodiscard]] int64_t score_sum() const{
        return qry_score + ref_score;
    }
    bool operator <(const PafDistance& rht) const{ // better
        if(*this == max()) return false;
        if(rht == max()) return true;
        assert(calc_sum_chk == rht.calc_sum_chk);
        if(cmp_mode == PafDistanceCompareMode::CALC_SUM_MODE) {
            if (score_sum() != rht.score_sum())
                return score_sum() < rht.score_sum();
        }else if(cmp_mode == PafDistanceCompareMode::QRY_SCORE_MODE){
            if(qry_score != rht.qry_score) return qry_score < rht.qry_score;
            if(ref_score != rht.ref_score) return ref_score < rht.ref_score;
        }else{
            assert(0); // never happens
        }
        if(anom != rht.anom) return anom < rht.anom;
        int64_t tot = qul_total ? qul_total: 1;
        int64_t rht_tot = rht.qul_total ? rht.qul_total : 1;
        return qul_nonzero * rht_tot > rht.qul_nonzero * tot;
    }
    bool operator >(const PafDistance& rht) const{
        return rht < *this;
    }
    bool operator ==(const PafDistance& rht) const{
        int64_t tot = qul_total ? qul_total: 1;
        int64_t rht_tot = rht.qul_total ? rht.qul_total : 1;
        return (qry_score == rht.qry_score) and (ref_score == rht.ref_score)
               and (anom == rht.anom) and (qul_nonzero * rht_tot == rht.qul_nonzero * tot);
    }
    bool operator <= (const PafDistance& rht) const{
        return *this == rht or *this < rht;
    }
    bool operator >= (const PafDistance &rht) const{
        return *this == rht or *this > rht;
    }
    bool operator !=(const PafDistance& rht) const{
        return not (*this == rht);
    }
    PafDistance operator+(const PafDistance& rht) const{
        assert(*this != max() and rht != max());
        assert(calc_sum_chk == rht.calc_sum_chk);
        return PafDistance(calc_sum_chk, qry_score + rht.qry_score,
                            ref_score + rht.ref_score, anom + rht.anom, qul_nonzero + rht.qul_nonzero, qul_total + rht.qul_total);
    }
    PafDistance operator-(const PafDistance& rht) const{
        assert(*this != max() and rht != max());
        assert(calc_sum_chk == rht.calc_sum_chk);
        return PafDistance(calc_sum_chk, qry_score - rht.qry_score, ref_score - rht.ref_score, anom - rht.anom, qul_nonzero - rht.qul_nonzero, qul_total - rht.qul_total);
    }
};

void get_overlap_range(PafReadData &paf_read_data, std::string_view cs_str);
PafEditData get_edited_paf_data(PafOutputData &paf_out, PafReadData &paf_read_data);
void solve_ctg_read(std::vector<PafReadData> &paf_ctg_data, std::vector<PafOutputData> &paf_ctg_out, std::vector<PafOutputData> &paf_ctg_alt_out, std::vector<std::vector<PafOutputData>> &paf_ctg_max_out);


#endif //ALIGNASM_PAF_DATA_HPP
