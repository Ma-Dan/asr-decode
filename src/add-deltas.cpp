#include "add-deltas.h"

DeltaFeatures::DeltaFeatures(const DeltaFeaturesOptions &opts): opts_(opts) {
  scales_.resize(opts.order+1);
  scales_[0].resize(1);
  scales_[0][0] = 1.0;  // trivial window for 0th order delta [i.e. baseline feats]

  for (int32 i = 1; i <= opts.order; i++) {
    vector<BaseFloat> &prev_scales = scales_[i-1],
        &cur_scales = scales_[i];
    int32 window = opts.window;  // this code is designed to still
    // work if instead we later make it an array and do opts.window[i-1],
    // or something like that. "window" is a parameter specifying delta-window
    // width which is actually 2*window + 1.
    int32 prev_offset = (static_cast<int32>(prev_scales.size()-1))/2,
        cur_offset = prev_offset + window;
    cur_scales.resize(prev_scales.size() + 2*window);  // also zeros it.

    BaseFloat normalizer = 0.0;
    for (int32 j = -window; j <= window; j++) {
      normalizer += j*j;
      for (int32 k = -prev_offset; k <= prev_offset; k++) {
        cur_scales[j+k+cur_offset] +=
            static_cast<BaseFloat>(j) * prev_scales[k+prev_offset];
      }
    }
    for(int32 i=0; i<cur_scales.size(); i++)
    {
        cur_scales[i] *= 1.0 / normalizer;
    }
  }
}

void DeltaFeatures::Process(const P_Matrix input_feats,
                            int32 frame,
                            BaseFloat *output_frame) const {
  int32 num_frames = input_feats->rows,
      feat_dim = input_feats->cols;
  for(int32 i=0; i<(opts_.order+1)*feat_dim; i++)
  {
      output_frame[i] = 0.0f;
  }
  for (int32 i = 0; i <= opts_.order; i++) {
    const vector<BaseFloat> &scales = scales_[i];
    int32 max_offset = (scales.size() - 1) / 2;
    BaseFloat* output = output_frame + i*feat_dim;
    for (int32 j = -max_offset; j <= max_offset; j++) {
      // if asked to read
      int32 offset_frame = frame + j;
      if (offset_frame < 0) offset_frame = 0;
      else if (offset_frame >= num_frames)
        offset_frame = num_frames - 1;
      BaseFloat scale = scales[j + max_offset];
      if (scale != 0.0)
      {
          for(int32 k=0; k<feat_dim; k++)
          {
              output[k] += scale * input_feats->data[offset_frame*input_feats->cols+k];
          }
      }
    }
  }
}

void ComputeDeltas(const DeltaFeaturesOptions &delta_opts,
                   const P_Matrix input_features,
                   P_Matrix output_features) {
  output_features->rows = input_features->rows;
  output_features->cols = input_features->cols*(delta_opts.order + 1);
  output_features->data.resize(output_features->rows * output_features->cols);
  DeltaFeatures delta(delta_opts);
  for (int32 r = 0; r < static_cast<int32>(input_features->rows); r++) {
    BaseFloat* row = output_features->data.data() + r*output_features->cols;
    delta.Process(input_features, r, row);
  }
}
