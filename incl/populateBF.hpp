
#pragma once


#include <data_structs.hpp>

SBF *populateBF(SEQLIST *assembly, string seed1, ulong aNum, double fpr, uint kMerSize);
ulong countAssemblyKmers(SEQLIST *assembly, uint kMerSize);
