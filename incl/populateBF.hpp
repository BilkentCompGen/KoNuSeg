#pragma once

#include <data_structs.hpp>

std::unique_ptr<SBF> populateBF(
    const SEQLIST& assembly,
    const std::string& seed1,
    uint64_t aNum,
    double fpr,
    uint32_t kMerSize,
    uint64_t& k_assembly_out);

uint64_t countAssemblyKmers(const SEQLIST& assembly, uint32_t kMerSize);
