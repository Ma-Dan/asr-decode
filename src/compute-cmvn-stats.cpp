#include "compute-cmvn-stats.h"

void InitCmvnStats(int32 dim, P_MatrixDouble stats) {
    stats->rows = 2;
    stats->cols = dim+1;
    stats->data.resize(2*(dim+1));
}

void AccCmvnStats(const BaseFloat* feats, int32 dim, BaseFloat weight, P_MatrixDouble stats) {
  // Remove these __restrict__ modifiers if they cause compilation problems.
  // It's just an optimization.
   double *__restrict__ mean_ptr = stats->data.data(),
       *__restrict__ var_ptr = stats->data.data()+stats->cols,
       *__restrict__ count_ptr = mean_ptr + dim;
   const BaseFloat * __restrict__ feats_ptr = feats;
  *count_ptr += weight;
  // Careful-- if we change the format of the matrix, the "mean_ptr < count_ptr"
  // statement below might become wrong.
  for (; mean_ptr < count_ptr; mean_ptr++, var_ptr++, feats_ptr++) {
    *mean_ptr += *feats_ptr * weight;
    *var_ptr +=  *feats_ptr * *feats_ptr * weight;
  }
}

void AccCmvnStats(const P_Matrix feats, P_MatrixDouble stats) {
  int32 num_frames = feats->rows;
  for (int32 i = 0; i < num_frames; i++) {
    const BaseFloat* this_frame = feats->data.data() + i * feats->cols;
    BaseFloat weight = 1.0;
    if (weight != 0.0)
      AccCmvnStats(this_frame, feats->cols, weight, stats);
  }
}

void ApplyCmvn(const P_MatrixDouble stats,
               bool var_norm,
               P_Matrix feats) {
  int32 dim = stats->cols - 1;

  double count = stats->data[dim];

  if (!var_norm) {
    vector<BaseFloat> offset;
    offset.resize(dim);
    for(int32 i=0; i<dim; i++)
    {
        offset[i] = -stats->data[i] / stats->data[dim];
    }
    for(int32 i=0; i<feats->rows; i++)
    {
        for(int32 j=0; j<feats->cols; j++)
        {
            feats->data[i*feats->cols+j] += offset[j];
        }
    }
    return;
  }
  // norm(0, d) = mean offset;
  // norm(1, d) = scale, e.g. x(d) <-- x(d)*norm(1, d) + norm(0, d).
  Matrix norm;
  norm.rows = 2;
  norm.cols = dim;
  norm.data.resize(2*dim);
  for (int32 d = 0; d < dim; d++) {
    double mean, offset, scale;
    mean = stats->data[d]/count;
    double var = (stats->data[1*stats->cols + d]/count) - mean*mean,
        floor = 1.0e-20;
    scale = 1.0 / sqrt(var);
    offset = -(mean*scale);
    norm.data[d] = offset;
    norm.data[1*norm.cols+d] = scale;
  }
  // Apply the normalization.
  //feats->MulColsVec(norm.Row(1));
  //feats->AddVecToRows(1.0, norm.Row(0));
}
