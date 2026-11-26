#pragma once

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <omp.h>

#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/data_stream.hpp>
#include <btllib/nthash.hpp>

namespace konuseg {

inline void trimRight(std::string& s) {
    const auto pos = s.find_last_not_of(" \t\n\r\f\v");
    if (pos == std::string::npos) s.clear();
    else s.erase(pos + 1);
}

inline void trim(std::string& s) {
    trimRight(s);
    const auto pos = s.find_first_not_of(" \t\n\r\f\v");
    if (pos != std::string::npos) s.erase(0, pos);
}

inline bool readLineFromStream(FILE* fp, std::string& line) {
    line.clear();
    if (fp == nullptr) return false;
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), fp) != nullptr) {
        const size_t len = std::strlen(buf);
        if (len > 0 && buf[len - 1] == '\n') {
            line.append(buf, len - 1);
            return true;
        }
        line.append(buf, len);
    }
    return !line.empty();
}

inline std::string deriveMergedBedPath(const std::string& outfile_path) {
    if (outfile_path.size() > 4 && outfile_path.compare(outfile_path.size() - 4, 4, ".bed") == 0) {
        return outfile_path.substr(0, outfile_path.size() - 4) + "_merged.bed";
    }
    return outfile_path + "_merged.bed";
}

} // namespace konuseg

class Config {
public:
    std::string infile_fastq;
    std::string infile_fasta;
    std::string outfile_bed;
    std::string save_bf_path;
    std::string load_bf_path;
    std::string save_cbf_path;
    std::string load_cbf_path;

    // Default values (used if the user doesn't provide them in the terminal)
    long long genome_size = 3000000000LL;
    std::string seed1 = "10010111001110100101110011101001";
    uint32_t k_len = 32;
    uint64_t k_num = 1500000000ULL;
    uint64_t a_num = 1500000000ULL;
    double fpr = 0.01;
    uint32_t window_size = 5;
    uint32_t seg_len = 500;
    uint32_t count_threshold = 1;

    void printConfig() const {
        std::cout << "\t\t[i]: Config file" << std::endl;
        std::cout << "\t\t[i]:" << infile_fastq << std::endl;
        std::cout << "\t\t[i]:" << infile_fasta << std::endl;
        std::cout << "\t\t[i]:" << outfile_bed << std::endl;
        std::cout << "\t\t[i]:" << genome_size << std::endl;
        std::cout << "\t\t[i]:" << seed1 << std::endl;
        std::cout << "\t\t[i]:" << k_len << std::endl;
        std::cout << "\t\t[i]:" << k_num << std::endl;
        std::cout << "\t\t[i]:" << a_num << std::endl;
        std::cout << "\t\t[i]:" << fpr << std::endl;
        std::cout << "\t\t[i]:" << window_size << std::endl;
        std::cout << "\t\t[i]:" << seg_len << std::endl;
        std::cout << "\t\t[i]:" << count_threshold << std::endl;
        if (!load_bf_path.empty()) std::cout << "\t\t[i]: load_bf=" << load_bf_path << std::endl;
        if (!save_bf_path.empty()) std::cout << "\t\t[i]: save_bf=" << save_bf_path << std::endl;
        if (!load_cbf_path.empty()) std::cout << "\t\t[i]: load_cbf=" << load_cbf_path << std::endl;
        if (!save_cbf_path.empty()) std::cout << "\t\t[i]: save_cbf=" << save_cbf_path << std::endl;
        std::cout << "\t\t[i]: //////////////////////" << std::endl;
    }
};

class Sequence {
public:
    std::string name;
    std::string data;

    Sequence() = default;
    Sequence(std::string n, std::string d)
        : name(std::move(n)), data(std::move(d)) {}
};

class ContigRepeatInfo {
public:
    std::string contigName;
    uint32_t start = 0;
    uint32_t end = 0;
    uint32_t numRepeats = 0;

    ContigRepeatInfo() = default;
    ContigRepeatInfo(std::string n, uint32_t s, uint32_t e, uint32_t r)
        : contigName(std::move(n)), start(s), end(e), numRepeats(r) {}

    void print() const {
        std::cout << contigName << "= s" << start << "-e" << end << ":" << numRepeats << std::endl;
    }
};

class SegmentInfo {
public:
    std::string segmentName;
    uint32_t start = 0;
    uint32_t end = 0;
    float count = 0.0f;

    SegmentInfo() = default;
    SegmentInfo(std::string n, uint32_t s, uint32_t e, float c)
        : segmentName(std::move(n)), start(s), end(e), count(c) {}
};

// Copy-number rounding heuristic used in the manuscript
inline float roundCount(float count) {
    const float intPart = std::floor(count);
    const float diff = count - intPart;

    if (count < 0.25f) {
        return 0.0f;
    } else if (count < 0.65f) {
        return 0.5f;
    } else {
        return (diff <= 0.64f) ? intPart : (intPart + 1.0f);
    }
}

using SEQLIST = std::vector<Sequence>;
using SBF = btllib::SeedBloomFilter;
using CBF = btllib::CountingBloomFilter16;
using CONTIGLIST = std::vector<ContigRepeatInfo>;
using SEGMENTINFO = std::vector<SegmentInfo>;
using SEGMENTTABLE = std::vector<SEGMENTINFO>;

struct PipelineState {
    uint64_t k_assembly = 0;
    uint64_t k_reads = 0;
    SEGMENTTABLE seg_table;
};
