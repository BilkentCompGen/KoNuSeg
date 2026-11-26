
#pragma once


#include <data_structs.hpp>


CBF *loadReads(
	string readsFile, SBF *seed_bf, ulong k_num, double fpr, string seed1, uint kMerSize
);

void populateCBF(
	CBF *cbf, const string& read, SBF *seed_bf, const string& seed1, uint kMerSize
);

void populate(
    const std::string& seq, const std::vector<std::string>& seeds, btllib::CountingBloomFilter16 *cbf
);