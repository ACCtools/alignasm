// /home/pentagon03/alignasm/HCC1954.p.paf -a /home/pentagon03/alignasm/alt.paf
#include "paf_data.hpp"

#include "csv.hpp"
#include "argparse/argparse.hpp"
#include <ankerl/unordered_dense.h>
#include <filesystem>

#include <iostream>
#include <cinttypes>
#include <string>
#include <utility>
#include <vector>
#include <cstring>
#include <cassert>
#include <charconv>
#include <string_view>

#ifdef NDEBUG
#include <tbb/parallel_for.h>
#endif


int32_t main(int argc, char** argv) {
    /** Test Session */
    argparse::ArgumentParser program("alignasm", "0.1.0");

    program.add_argument("PAF_LOC")
            .help("Location of PAF file")
            .required()
            .nargs(1);

    program.add_argument("-t", "--thread")
            .help("Number of threads")
            .default_value(1)
            .scan<'d', int>()
            .metavar("THREAD");


    program.add_argument("-a", "--alt")
            .help("Location of alternative PAF file")
            .nargs(1)
            .metavar("PAF_ALT_LOC");

    program.add_argument("-b", "--alt_baseline")
            .help("Baseline for coverage of alternative PAF file")
            .default_value(0.5)
            .scan<'f', double>()
            .metavar("ALT_BASELINE");

    try {
        program.parse_args(argc, argv);
    }
    catch (...) {
        std::cerr << program;
        return 1;
    }

    std::filesystem::path paf_loc {program.get<std::string>("PAF_LOC")};
#ifdef NDEBUG
    int num_thread = program.get<int>("--thread");
#endif
    if (paf_loc.extension() != ".paf") {
        std::cerr << "Wrong PAF file : " << std::filesystem::absolute(paf_loc);
        std::cerr << program;
        return 1;
    }

    /** CSV read */
    csv::CSVFormat format;
    format.delimiter('\t')
            .variable_columns(true)
            .variable_columns(csv::VariableColumnPolicy::KEEP)
            .no_header();

    std::string filename = std::filesystem::absolute(paf_loc);
    csv::CSVReader reader(filename, format);

    /** CSV Process Data */
    using ChrMap = ankerl::unordered_dense::map <std::string, int32_t>;
    using ChrRevMap = ankerl::unordered_dense::map <int32_t, std::string >;

    ChrMap chr_map{};
    ChrRevMap chr_rev_map {};
    ChrMap paf_map {};
    int32_t chr_map_index = 0;

    std::vector<std::vector<PafReadData>> paf_data;
    std::vector<PafReadData> ctg_data_vector;
    std::vector<std::string> ctg_name_vector;
    std::string ctg_chr;
    for (int32_t ctg_index = 0, paf_index = 0; csv::CSVRow& read : reader) {
        std::string qry_chr, ref_chr;
        ref_chr = read[PAF_REF_CHR].get<std::string>();
        qry_chr = read[PAF_QRY_CHR].get<std::string>();

        if (ctg_chr.empty()) {
            ctg_chr = qry_chr;
        }

        if (not chr_map.contains(ref_chr)) {
            chr_map[ref_chr] = chr_map_index;
            chr_rev_map[chr_map_index] = ref_chr;
            chr_map_index++;
        }

        if (ctg_chr != qry_chr) {
            paf_data.push_back(ctg_data_vector);
            ctg_name_vector.push_back(ctg_chr);
            ctg_chr = qry_chr;
            ctg_index = 0;
            ctg_data_vector.clear();

            paf_index += 1;
        }

        PafReadData paf_read_data{};
        paf_map[qry_chr] = paf_index;
        paf_read_data.paf_index = paf_index;
        paf_read_data.ctg_index = ctg_index;
        ctg_index += 1; // starts from 0

        paf_read_data.qry_total_length = read[PAF_QRY_TOT].get<int64_t>();
        paf_read_data.qry_str = read[PAF_QRY_STR].get<int64_t>();
        paf_read_data.qry_end = read[PAF_QRY_END].get<int64_t>();
        paf_read_data.qry_end--;
        assert(paf_read_data.qry_str <= paf_read_data.qry_end); // start <= end

        paf_read_data.ref_total_length = read[PAF_REF_TOT].get<int64_t>();
        paf_read_data.ref_str = read[PAF_REF_STR].get<int64_t>();
        paf_read_data.ref_end = read[PAF_REF_END].get<int64_t>();
        paf_read_data.ref_end--;
        assert(paf_read_data.ref_str <= paf_read_data.ref_end); // start <= end

        paf_read_data.ref_chr = chr_map[ref_chr];

        paf_read_data.aln_fwd = read[PAF_ALN_FWD].get<std::string>()[0] == '+';
        if(not paf_read_data.aln_fwd) {
            std::swap(paf_read_data.ref_str, paf_read_data.ref_end);
            assert(paf_read_data.ref_str >= paf_read_data.ref_end); // start >= end
        }

        paf_read_data.map_qul = read[PAF_MAT_QUL].get<uint8_t>();


        // cs tag is last in alignmentpaf
        get_overlap_range(paf_read_data, read[read.size() - 1].get<std::string_view>());
        ctg_data_vector.push_back(paf_read_data);
    }

    ctg_name_vector.push_back(ctg_chr);
    paf_data.push_back(ctg_data_vector);

    assert(not ctg_data_vector.empty() and "data should not be empty");


    if (program.is_used("--alt")) {
        double ALT_BASELINE = program.get<double>("--alt_baseline");

        std::filesystem::path paf_alt_loc {program.get<std::string>("--alt")};
        if (paf_alt_loc.extension() != ".paf") {
            std::cerr << "Wrong PAF file : " << std::filesystem::absolute(paf_loc);
            std::cerr << program;
            return 1;
        }

        filename = std::filesystem::absolute(paf_alt_loc);
        csv::CSVReader alt_reader(filename, format);

        auto parseString = [](const std::string& input) -> std::pair<std::string, int64_t> {
            // ':'를 기준으로 문자열을 분리
            size_t pos = input.find(':');
            if (pos == std::string::npos) {
                throw std::invalid_argument("Invalid input string format");
            }

            // 앞부분을 추출 (ptg000002l)
            std::string firstPart = input.substr(0, pos);

            // ':' 다음 부분을 추출하고 '-' 이전까지의 숫자를 추출
            size_t start = pos + 1;
            size_t end = input.find('-', start);
            if (end == std::string::npos) {
                end = input.length();
            }

            int64_t secondPart;
            auto result = std::from_chars(input.data() + start, input.data() + end, secondPart);
            if (result.ec != std::errc()) {
                throw std::invalid_argument("Error parsing number");
            }

            return std::make_pair(firstPart, secondPart - 1);
        };


        std::string tar_real_qry_chr = {};
        int64_t tar_qry_offset = -1;
        bool tar_flag = false;
        double tar_ratio = 0;
        double aln_ratio;
        PafReadData ratio_max_paf_data{};

        ctg_chr.clear();
        for (csv::CSVRow& read : alt_reader) {
            std::string qry_chr, ref_chr;
            ref_chr = read[PAF_REF_CHR].get<std::string>();
            qry_chr = read[PAF_QRY_CHR].get<std::string>();

            if (not chr_map.contains(ref_chr)) {
                chr_map[ref_chr] = chr_map_index;
                chr_rev_map[chr_map_index] = ref_chr;
                chr_map_index++;
            }

            auto [real_qry_chr, qry_offset] = parseString(qry_chr);
            auto& ctg_last_data = paf_data[paf_map[real_qry_chr]].back();

            PafReadData paf_read_data{};
            paf_read_data.paf_index = ctg_last_data.paf_index;
            paf_read_data.ctg_index = ctg_last_data.ctg_index + 1;

            paf_read_data.qry_total_length = ctg_last_data.qry_total_length;
            paf_read_data.qry_str = read[PAF_QRY_STR].get<int64_t>() + qry_offset;
            paf_read_data.qry_end = read[PAF_QRY_END].get<int64_t>() + qry_offset;
            paf_read_data.qry_end--;
            assert(paf_read_data.qry_str <= paf_read_data.qry_end); // start <= end

            paf_read_data.ref_total_length = read[PAF_REF_TOT].get<int64_t>();
            paf_read_data.ref_str = read[PAF_REF_STR].get<int64_t>();
            paf_read_data.ref_end = read[PAF_REF_END].get<int64_t>();
            paf_read_data.ref_end--;
            assert(paf_read_data.ref_str <= paf_read_data.ref_end); // start <= end

            paf_read_data.ref_chr = chr_map[ref_chr];

            paf_read_data.aln_fwd = read[PAF_ALN_FWD].get<std::string>()[0] == '+';
            if(not paf_read_data.aln_fwd) {
                std::swap(paf_read_data.ref_str, paf_read_data.ref_end);
                assert(paf_read_data.ref_str >= paf_read_data.ref_end); // start >= end
            }

            paf_read_data.map_qul = read[PAF_MAT_QUL].get<uint8_t>();

            // cs tag is last in alignment
            get_overlap_range(paf_read_data, read[read.size() - 1].get<std::string_view>());

            if (tar_qry_offset == -1) {
                tar_qry_offset = qry_offset;
                tar_real_qry_chr = real_qry_chr;
            }

            if (tar_qry_offset != qry_offset or tar_real_qry_chr != real_qry_chr) {
                if (not tar_flag) {
                    paf_data[paf_map[tar_real_qry_chr]].push_back(ratio_max_paf_data);
                }

                tar_flag = false;
                tar_ratio = 0;
                tar_qry_offset = qry_offset;
                tar_real_qry_chr = real_qry_chr;
            }

            aln_ratio = read[PAF_ALN_LEN].get<double>() / read[PAF_QRY_TOT].get<double>();

            if (aln_ratio > tar_ratio) {
                tar_ratio = aln_ratio;
                ratio_max_paf_data = paf_read_data;
            }

            if (aln_ratio > ALT_BASELINE) {
                paf_data[paf_map[real_qry_chr]].push_back(paf_read_data);
                tar_flag = true;
            }
        }
    }

    std::cout << "File read complete!\n";

    /** Output Data */
    std::vector<std::vector<PafOutputData>> paf_out_data(paf_data.size()), paf_alt_out_data(paf_data.size());
    std::vector< std::vector<std::vector<PafOutputData>> > paf_max_out_datas(paf_data.size());

#ifdef NDEBUG
    if (num_thread > 1) {
        std::cout << "Task in Multi-Threading" << '\n';
        tbb::task_arena arena(num_thread);

        arena.execute([&paf_data, &paf_out_data, &paf_alt_out_data, &paf_max_out_datas] {
            tbb::parallel_for(tbb::blocked_range<unsigned long>(0, paf_data.size()),
                              [&paf_data, &paf_out_data, &paf_alt_out_data, &paf_max_out_datas](const tbb::blocked_range<unsigned long>& range) {
                                  for (auto i = range.begin(); i < range.end(); i++) {
                                      solve_ctg_read(paf_data[i], paf_out_data[i], paf_alt_out_data[i], paf_max_out_datas[i]);
                                  }
                              });
        });
    } else {
        for (int32_t i = 0; i < paf_data.size(); i++) {
            std::cout << "Task " << (i+1) << " in Progress." << '\n';
            solve_ctg_read(paf_data[i], paf_out_data[i], paf_alt_out_data[i], paf_max_out_datas[i]);
        }
    }
#else
    for (int32_t i = 0; i < paf_data.size(); i++) {
        solve_ctg_read(paf_data[i], paf_out_data[i], paf_alt_out_data[i], paf_max_out_datas[i]);
    }
#endif

    auto process_output = [&](auto &paf_out_data_, const std::string &prefix = "") {
        /* Write Output */
        std::ofstream paf_out(std::filesystem::absolute(paf_loc).replace_extension(".aln" + prefix + ".paf"));
        auto writer = csv::make_tsv_writer(paf_out);

        assert(ctg_name_vector.size() == paf_out_data_.size() and paf_out_data_.size() == paf_data.size());

        auto ctg_n = ctg_name_vector.size();
        for (auto i = 0; i < ctg_n; i++) {
            for (auto& paf_line : paf_out_data_[i]) {
                auto& paf_data_line = paf_data[i][paf_line.ctg_index];

                writer << std::vector<std::string>({ctg_name_vector[i], std::to_string(paf_data_line.qry_total_length), std::to_string(paf_line.edited_qry_str), std::to_string(paf_line.edited_qry_end + 1),
                                                    paf_data_line.aln_fwd ? "+" : "-", chr_rev_map[paf_data_line.ref_chr],
                                                    std::to_string(paf_data_line.ref_total_length), std::to_string(paf_line.edited_ref_str), std::to_string(paf_line.edited_ref_end + 1)});
            }
        }
    };

    std::cout << "Writing Output" << '\n';
    process_output(paf_out_data);
    process_output(paf_alt_out_data, ".alt");

}
