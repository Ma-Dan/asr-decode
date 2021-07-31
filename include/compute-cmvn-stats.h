#ifndef COMPUTE_CMVN_STATS
#define COMPUTE_CMVN_STATS

#include "common.h"

void InitCmvnStats(int32 dim, P_MatrixDouble stats);
void AccCmvnStats(const P_Matrix feats, P_MatrixDouble stats);
void ApplyCmvn(const P_MatrixDouble stats, bool var_norm, P_Matrix feats);

#endif