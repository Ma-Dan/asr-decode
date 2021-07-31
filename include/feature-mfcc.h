#ifndef FEATURE_MFCC
#define FEATURE_MFCC

#include "common.h"
#include "srfft.h"

int32 RoundUpToNearestPowerOfTwo(int32 n);

enum WindowsType
{
    hanning = 0,
    sine,
    hamming,
    povey,
    rectangular,
    blackman
};

struct FrameExtractionOptions
{
    BaseFloat samp_freq;
    BaseFloat frame_shift_ms;  // in milliseconds.
    BaseFloat frame_length_ms;  // in milliseconds.
    BaseFloat dither;  // Amount of dithering, 0.0 means no dither.
    BaseFloat preemph_coeff;  // Preemphasis coefficient.
    bool remove_dc_offset;  // Subtract mean of wave before FFT.
    vector<BaseFloat> window;
    WindowsType window_type;
    BaseFloat blackman_coeff;

    FrameExtractionOptions():
      samp_freq(8000),
      frame_shift_ms(10.0),
      frame_length_ms(25.0),
      dither(1.0),
      preemph_coeff(0.97),
      remove_dc_offset(true),
      window_type(povey),
      blackman_coeff(0.42){};

    int32 WindowShift() const {
        return static_cast<int32>(samp_freq * 0.001 * frame_shift_ms);
    }
    int32 WindowSize() const {
        return static_cast<int32>(samp_freq * 0.001 * frame_length_ms);
    }
    int32 PaddedWindowSize() const {
        return RoundUpToNearestPowerOfTwo(WindowSize());
    }
};

struct MelBanksOptions {
  int32 num_bins;  // e.g. 25; number of triangular bins
  BaseFloat low_freq;  // e.g. 20; lower frequency cutoff
  BaseFloat high_freq;  // an upper frequency cutoff; 0 -> no cutoff, negative
  // ->added to the Nyquist frequency to get the cutoff.
  BaseFloat vtln_low;  // vtln lower cutoff of warping function.
  BaseFloat vtln_high;  // vtln upper cutoff of warping function: if negative, added
                        // to the Nyquist frequency to get the cutoff.

  explicit MelBanksOptions(int num_bins = 25)
      : num_bins(num_bins), low_freq(20), high_freq(0), vtln_low(100),
        vtln_high(-500) {}
};

class MelBanks {
 public:

  static inline BaseFloat InverseMelScale(BaseFloat mel_freq) {
    return 700.0f * (expf (mel_freq / 1127.0f) - 1.0f);
  }

  static inline BaseFloat MelScale(BaseFloat freq) {
    return 1127.0f * logf (1.0f + freq / 700.0f);
  }

  static BaseFloat VtlnWarpFreq(BaseFloat vtln_low_cutoff,
                                BaseFloat vtln_high_cutoff,  // discontinuities in warp func
                                BaseFloat low_freq,
                                BaseFloat high_freq,  // upper+lower frequency cutoffs in
                                // the mel computation
                                BaseFloat vtln_warp_factor,
                                BaseFloat freq);

  static BaseFloat VtlnWarpMelFreq(BaseFloat vtln_low_cutoff,
                                   BaseFloat vtln_high_cutoff,
                                   BaseFloat low_freq,
                                   BaseFloat high_freq,
                                   BaseFloat vtln_warp_factor,
                                   BaseFloat mel_freq);


  MelBanks(const MelBanksOptions &opts,
           const FrameExtractionOptions &frame_opts,
           BaseFloat vtln_warp_factor);

  /// Compute Mel energies (note: not log enerties).
  /// At input, "fft_energies" contains the FFT energies (not log).
  void Compute(const vector<BaseFloat> &fft_energies,
               vector<BaseFloat> &mel_energies_out) const;

  int32 NumBins() const { return bins_.size(); }

  // returns vector of central freq of each bin; needed by plp code.
  const vector<BaseFloat> &GetCenterFreqs() const { return center_freqs_; }

  const std::vector<std::pair<int32, vector<BaseFloat> > >& GetBins() const {
    return bins_;
  }

 private:

  // center frequencies of bins, numbered from 0 ... num_bins-1.
  // Needed by GetCenterFreqs().
  vector<BaseFloat> center_freqs_;

  // the "bins_" vector is a vector, one for each bin, of a pair:
  // (the first nonzero fft-bin), (the vector of weights).
  std::vector<std::pair<int32, vector<BaseFloat> > > bins_;
};

struct MfccOptions
{
    MelBanksOptions mel_opts;
    int num_ceps;   // e.g. 13: num cepstral coeffs, counting zero.
    BaseFloat cepstral_lifter;  // Scaling factor on cepstra for HTK compatibility.
                              // if 0.0, no liftering is done.

    MfccOptions()
    : mel_opts(23),
      num_ceps(13),
      cepstral_lifter(22.0)
    {};
};

struct RandomState {
    RandomState();
    unsigned seed;
};

class MfccComputer
{
    public:
        MfccComputer();
        ~MfccComputer();
        void ComputeFeatures(const vector<BaseFloat> &wave, BaseFloat sample_freq, BaseFloat vtln_warp, vector<BaseFloat* > &output);

    private:
        int32 NumFrames(int64 num_samples);
        void ExtractWindow(const vector<BaseFloat> &wave, int32 f, BaseFloat vtln_warp, vector<BaseFloat> &window, BaseFloat* output);
        void ProcessWindow(vector<BaseFloat> window, BaseFloat vtln_warp, BaseFloat* output);
        const MelBanks *GetMelBanks(BaseFloat vtln_warp);

        MfccOptions mfccOptions;
        FrameExtractionOptions frameOptions;
        std::map<BaseFloat, MelBanks*> mel_banks_;  // BaseFloat is VTLN coefficient.
        SplitRadixRealFft *srfft;

        vector<BaseFloat> mel_energies_;
        Matrix dct_matrix;
        vector<BaseFloat> lifter_coeffs_;
};

#endif