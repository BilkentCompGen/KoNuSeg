#pragma once

#include <data_structs.hpp>

bool checkRepeats(
    const SEQLIST& assembly,
    const CBF* cbf,
    uint32_t kMerSize,
    const std::string& seed1,
    uint32_t windowSize,
    const std::string& outfile_bed,
    uint32_t seg_len,
    uint32_t count_threshold,
    PipelineState& state);
