#pragma once

#include <data_structs.hpp>

bool mergeSegments(
    const std::string& outfile_bed,
    uint32_t k_len,
    uint32_t seg_len,
    uint32_t count_threshold,
    const SEGMENTTABLE& seg_table);
