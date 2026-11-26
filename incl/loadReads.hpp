#pragma once

#include <data_structs.hpp>

std::unique_ptr<CBF> loadReads(
    const std::string& readsFile,
    const SBF* seed_bf,
    uint64_t k_num,
    double fpr,
    const std::string& seed1,
    uint32_t kMerSize,
    uint64_t k_assembly,
    uint64_t& k_reads_out);

uint64_t populateCBF(
    CBF* cbf,
    const std::string& read,
    const SBF* seed_bf,
    const std::string& seed1,
    uint32_t kMerSize);

void populate(
    const std::string& seq,
    const std::vector<std::string>& seeds,
    btllib::CountingBloomFilter16* cbf);
