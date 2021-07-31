#ifndef ADD_DELTAS
#define ADD_DELTAS

#include "common.h"

struct DeltaFeaturesOptions {
  int32 order;
  int32 window;  // e.g. 2; controls window size (window size is 2*window + 1)
  // the behavior at the edges is to replicate the first or last frame.
  // this is not configurable.

  DeltaFeaturesOptions(int32 order = 2, int32 window = 2):
      order(order), window(window) { }
};

class DeltaFeatures {
 public:
  // This class provides a low-level function to compute delta features.
  // The function takes as input a matrix of features and a frame index
  // that it should compute the deltas on.  It puts its output in an object
  // of type VectorBase, of size (original-feature-dimension) * (opts.order+1).
  // This is not the most efficient way to do the computation, but it's
  // state-free and thus easier to understand

  explicit DeltaFeatures(const DeltaFeaturesOptions &opts);

  void Process(const P_Matrix input_feats,
               int32 frame,
               BaseFloat *output_frame) const;
 private:
  DeltaFeaturesOptions opts_;
  std::vector<vector<BaseFloat> > scales_;  // a scaling window for each
  // of the orders, including zero: multiply the features for each
  // dimension by this window.
};

void ComputeDeltas(const DeltaFeaturesOptions &delta_opts,
                   const P_Matrix input_features,
                   P_Matrix output_features);

#endif