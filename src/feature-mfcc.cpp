#include "feature-mfcc.h"

int32 RoundUpToNearestPowerOfTwo(int32 n) {
        n--;
        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;
        return n+1;
}

int Rand(struct RandomState* state) {
    if (state) {
        return rand_r(&(state->seed));
    } else {
        return rand();
    }
}

RandomState::RandomState() {
  // we initialize it as Rand() + 27437 instead of just Rand(), because on some
  // systems, e.g. at the very least Mac OSX Yosemite and later, it seems to be
  // the case that rand_r when initialized with rand() will give you the exact
  // same sequence of numbers that rand() will give if you keep calling rand()
  // after that initial call.  This can cause problems with repeated sequences.
  // For example if you initialize two RandomState structs one after the other
  // without calling rand() in between, they would give you the same sequence
  // offset by one (if we didn't have the "+ 27437" in the code).  27437 is just
  // a randomly chosen prime number.
  seed = unsigned(Rand(NULL)) + 27437;
}

/// Returns a random number strictly between 0 and 1.
inline float RandUniform(struct RandomState* state = NULL) {
  return static_cast<float>((Rand(state) + 1.0) / (RAND_MAX+2.0));
}

inline float RandGauss(struct RandomState* state = NULL) {
  return static_cast<float>(sqrtf (-2 * logf(RandUniform(state)))
                            * cosf(2*M_PI*RandUniform(state)));
}

void Dither(vector<BaseFloat> &waveform, int32 frame_length, BaseFloat dither_value) {
    if (dither_value == 0.0)
    {
        return;
    }
    BaseFloat *data = waveform.data();
    RandomState rstate;
    for (int32 i = 0; i < frame_length; i++)
    {
        data[i] += RandGauss(&rstate) * dither_value;
    }
}

BaseFloat Sum(vector<BaseFloat> window)
{
    BaseFloat sum = 0.0f;

    for(int i=0; i<window.size(); i++)
    {
        sum += window[i];
    }

    return sum;
}

void Preemphasize(vector<BaseFloat> &waveform, int32 frame_length, BaseFloat preemph_coeff)
{
    if (preemph_coeff == 0.0)
    {
        return;
    }
    for (int32 i = frame_length-1; i > 0; i--)
    {
        waveform[i] -= preemph_coeff * waveform[i-1];
    }

    waveform[0] -= preemph_coeff * waveform[0];
}

void ComputePowerSpectrum(vector<BaseFloat> &waveform) {
  int32 dim = waveform.size();

  // no, letting it be non-power-of-two for now.
  // KALDI_ASSERT(dim > 0 && (dim & (dim-1) == 0));  // make sure a power of two.. actually my FFT code
  // does not require this (dan) but this is better in case we use different code [dan].

  // RealFft(waveform, true);  // true == forward (not inverse) FFT; makes no difference here,
  // as we just want power spectrum.

  // now we have in waveform, first half of complex spectrum
  // it's stored as [real0, realN/2, real1, im1, real2, im2, ...]
  int32 half_dim = dim/2;
  BaseFloat first_energy = waveform[0] * waveform[0],
      last_energy = waveform[1] * waveform[1];  // handle this special case
  for (int32 i = 1; i < half_dim; i++) {
    BaseFloat real = waveform[i*2], im = waveform[i*2 + 1];
    waveform[i] = real*real + im*im;
  }
  waveform[0] = first_energy;
  waveform[half_dim] = last_energy;  // Will actually never be used, and anyway
  // if the signal has been bandlimited sensibly this should be zero.
}

void ApplyFloor(vector<BaseFloat> &v, BaseFloat floor_val)
{
    for (int32 i = 0; i < v.size(); i++) {
      v[i] = std::max(v[i], floor_val);
    }
}

void ApplyLog(vector<BaseFloat> &v)
{
    for (int32 i = 0; i < v.size(); i++) {
      v[i] = logf(v[i]);
    }
}

void PrepareMatrix(P_Matrix m, int32 rows, int32 cols)
{
    m->rows = rows;
    m->cols = cols;

    m->data.resize(rows * cols);
}

void ComputeDctMatrix(P_Matrix M) {
  //KALDI_ASSERT(M->NumRows() == M->NumCols());
  int32 K = M->rows;
  int32 N = M->cols;

  BaseFloat normalizer = sqrt(1.0 / static_cast<BaseFloat>(N));  // normalizer for
  // X_0.
  for (int32 j = 0; j < N; j++) M->data[0*M->cols + j] = normalizer;
  normalizer = sqrt(2.0 / static_cast<BaseFloat>(N));  // normalizer for other
   // elements.
  for (int32 k = 1; k < K; k++)
    for (int32 n = 0; n < N; n++)
      M->data[k*M->cols + n] = normalizer
          * cos( static_cast<double>(M_PI)/N * (n + 0.5) * k );
}

void ComputeLifterCoeffs(BaseFloat Q, vector<BaseFloat> &coeffs) {
  // Compute liftering coefficients (scaling on cepstral coeffs)
  // coeffs are numbered slightly differently from HTK: the zeroth
  // index is C0, which is not affected.
  for (int32 i = 0; i < coeffs.size(); i++)
    coeffs[i] = 1.0 + 0.5 * Q * sin (M_PI * i / Q);
}

void PrepareFeatureWindowFunction(FrameExtractionOptions &opts) {
  int32 frame_length = opts.WindowSize();
  opts.window.resize(frame_length);
  double a = M_2PI / (frame_length-1);
  for (int32 i = 0; i < frame_length; i++) {
    double i_fl = static_cast<double>(i);
    if (opts.window_type == hanning) {
      opts.window[i] = 0.5  - 0.5*cos(a * i_fl);
    } else if (opts.window_type == sine) {
      // when you are checking ws wikipedia, please
      // note that 0.5 * a = M_PI/(frame_length-1)
      opts.window[i] = sin(0.5 * a * i_fl);
    } else if (opts.window_type == hamming) {
      opts.window[i] = 0.54 - 0.46*cos(a * i_fl);
    } else if (opts.window_type == povey) {  // like hamming but goes to zero at edges.
      opts.window[i] = pow(0.5 - 0.5*cos(a * i_fl), 0.85);
    } else if (opts.window_type == rectangular) {
      opts.window[i] = 1.0;
    } else if (opts.window_type == blackman) {
      opts.window[i] = opts.blackman_coeff - 0.5*cos(a * i_fl) +
        (0.5 - opts.blackman_coeff) * cos(2 * a * i_fl);
    }
  }
}

MelBanks::MelBanks(const MelBanksOptions &opts,
                   const FrameExtractionOptions &frame_opts,
                   BaseFloat vtln_warp_factor) {
  int32 num_bins = opts.num_bins;
  BaseFloat sample_freq = frame_opts.samp_freq;
  int32 window_length_padded = frame_opts.PaddedWindowSize();
  int32 num_fft_bins = window_length_padded / 2;
  BaseFloat nyquist = 0.5 * sample_freq;

  BaseFloat low_freq = opts.low_freq, high_freq;
  if (opts.high_freq > 0.0)
    high_freq = opts.high_freq;
  else
    high_freq = nyquist + opts.high_freq;

  BaseFloat fft_bin_width = sample_freq / window_length_padded;
  // fft-bin width [think of it as Nyquist-freq / half-window-length]

  BaseFloat mel_low_freq = MelScale(low_freq);
  BaseFloat mel_high_freq = MelScale(high_freq);

  // divide by num_bins+1 in next line because of end-effects where the bins
  // spread out to the sides.
  BaseFloat mel_freq_delta = (mel_high_freq - mel_low_freq) / (num_bins+1);

  BaseFloat vtln_low = opts.vtln_low,
      vtln_high = opts.vtln_high;
  if (vtln_high < 0.0) {
    vtln_high += nyquist;
  }

  bins_.resize(num_bins);
  center_freqs_.resize(num_bins);

  for (int32 bin = 0; bin < num_bins; bin++) {
    BaseFloat left_mel = mel_low_freq + bin * mel_freq_delta,
        center_mel = mel_low_freq + (bin + 1) * mel_freq_delta,
        right_mel = mel_low_freq + (bin + 2) * mel_freq_delta;

    if (vtln_warp_factor != 1.0) {
      left_mel = VtlnWarpMelFreq(vtln_low, vtln_high, low_freq, high_freq,
                                 vtln_warp_factor, left_mel);
      center_mel = VtlnWarpMelFreq(vtln_low, vtln_high, low_freq, high_freq,
                                 vtln_warp_factor, center_mel);
      right_mel = VtlnWarpMelFreq(vtln_low, vtln_high, low_freq, high_freq,
                                  vtln_warp_factor, right_mel);
    }
    center_freqs_[bin] = InverseMelScale(center_mel);
    // this_bin will be a vector of coefficients that is only
    // nonzero where this mel bin is active.
    vector<BaseFloat> this_bin(num_fft_bins);
    int32 first_index = -1, last_index = -1;
    for (int32 i = 0; i < num_fft_bins; i++) {
      BaseFloat freq = (fft_bin_width * i);  // Center frequency of this fft
                                             // bin.
      BaseFloat mel = MelScale(freq);
      if (mel > left_mel && mel < right_mel) {
        BaseFloat weight;
        if (mel <= center_mel)
          weight = (mel - left_mel) / (center_mel - left_mel);
        else
         weight = (right_mel-mel) / (right_mel-center_mel);
        this_bin[i] = weight;
        if (first_index == -1)
          first_index = i;
        last_index = i;
      }
    }

    bins_[bin].first = first_index;
    int32 size = last_index + 1 - first_index;
    bins_[bin].second.resize(size);
    for(int32 i=0; i<size; i++)
    {
        bins_[bin].second[i] = this_bin[first_index+i];
    }
  }
}

BaseFloat MelBanks::VtlnWarpFreq(BaseFloat vtln_low_cutoff,  // upper+lower frequency cutoffs for VTLN.
                                 BaseFloat vtln_high_cutoff,
                                 BaseFloat low_freq,  // upper+lower frequency cutoffs in mel computation
                                 BaseFloat high_freq,
                                 BaseFloat vtln_warp_factor,
                                 BaseFloat freq) {
  /// This computes a VTLN warping function that is not the same as HTK's one,
  /// but has similar inputs (this function has the advantage of never producing
  /// empty bins).

  /// This function computes a warp function F(freq), defined between low_freq and
  /// high_freq inclusive, with the following properties:
  ///  F(low_freq) == low_freq
  ///  F(high_freq) == high_freq
  /// The function is continuous and piecewise linear with two inflection
  ///   points.
  /// The lower inflection point (measured in terms of the unwarped
  ///  frequency) is at frequency l, determined as described below.
  /// The higher inflection point is at a frequency h, determined as
  ///   described below.
  /// If l <= f <= h, then F(f) = f/vtln_warp_factor.
  /// If the higher inflection point (measured in terms of the unwarped
  ///   frequency) is at h, then max(h, F(h)) == vtln_high_cutoff.
  ///   Since (by the last point) F(h) == h/vtln_warp_factor, then
  ///   max(h, h/vtln_warp_factor) == vtln_high_cutoff, so
  ///   h = vtln_high_cutoff / max(1, 1/vtln_warp_factor).
  ///     = vtln_high_cutoff * min(1, vtln_warp_factor).
  /// If the lower inflection point (measured in terms of the unwarped
  ///   frequency) is at l, then min(l, F(l)) == vtln_low_cutoff
  ///   This implies that l = vtln_low_cutoff / min(1, 1/vtln_warp_factor)
  ///                       = vtln_low_cutoff * max(1, vtln_warp_factor)


  if (freq < low_freq || freq > high_freq) return freq;  // in case this gets called
  // for out-of-range frequencies, just return the freq.

  BaseFloat one = 1.0;
  BaseFloat l = vtln_low_cutoff * std::max(one, vtln_warp_factor);
  BaseFloat h = vtln_high_cutoff * std::min(one, vtln_warp_factor);
  BaseFloat scale = 1.0 / vtln_warp_factor;
  BaseFloat Fl = scale * l;  // F(l);
  BaseFloat Fh = scale * h;  // F(h);
  // slope of left part of the 3-piece linear function
  BaseFloat scale_left = (Fl - low_freq) / (l - low_freq);
  // [slope of center part is just "scale"]

  // slope of right part of the 3-piece linear function
  BaseFloat scale_right = (high_freq - Fh) / (high_freq - h);

  if (freq < l) {
    return low_freq + scale_left * (freq - low_freq);
  } else if (freq < h) {
    return scale * freq;
  } else {  // freq >= h
    return high_freq + scale_right * (freq - high_freq);
  }
}

BaseFloat MelBanks::VtlnWarpMelFreq(BaseFloat vtln_low_cutoff,  // upper+lower frequency cutoffs for VTLN.
                                    BaseFloat vtln_high_cutoff,
                                    BaseFloat low_freq,  // upper+lower frequency cutoffs in mel computation
                                    BaseFloat high_freq,
                                    BaseFloat vtln_warp_factor,
                                    BaseFloat mel_freq) {
  return MelScale(VtlnWarpFreq(vtln_low_cutoff, vtln_high_cutoff,
                               low_freq, high_freq,
                               vtln_warp_factor, InverseMelScale(mel_freq)));
}


// "power_spectrum" contains fft energies.
void MelBanks::Compute(const vector<BaseFloat> &power_spectrum,
                       vector<BaseFloat> &mel_energies_out) const {
  int32 num_bins = bins_.size();

  for (int32 i = 0; i < num_bins; i++) {
    int32 offset = bins_[i].first;
    const vector<BaseFloat> &v(bins_[i].second);
    BaseFloat energy = 0.0f;
    for(int32 j=0; j<v.size(); j++)
    {
        energy += v[j] * power_spectrum[offset+j];
    }
    mel_energies_out[i] = energy;
  }
}

MfccComputer::MfccComputer()
{
    PrepareFeatureWindowFunction(frameOptions);
    srfft = new SplitRadixRealFft(256);

    GetMelBanks(1.0);
    mel_energies_.resize(mfccOptions.mel_opts.num_bins);

    int32 num_bins = mfccOptions.mel_opts.num_bins;
    Matrix tmp_dct_matrix;
    PrepareMatrix(&tmp_dct_matrix, num_bins, num_bins);
    ComputeDctMatrix(&tmp_dct_matrix);
    PrepareMatrix(&dct_matrix, mfccOptions.num_ceps, num_bins);
    for(int32 i=0; i<dct_matrix.rows; i++)
    {
        for(int32 j=0; j<dct_matrix.cols; j++)
        {
            dct_matrix.data[i*dct_matrix.cols + j] = tmp_dct_matrix.data[i*tmp_dct_matrix.cols + j];
        }
    }

    if (mfccOptions.cepstral_lifter != 0.0) {
        lifter_coeffs_.resize(mfccOptions.num_ceps);
        ComputeLifterCoeffs(mfccOptions.cepstral_lifter, lifter_coeffs_);
    }
}

MfccComputer::~MfccComputer()
{
    delete srfft;
}

void MfccComputer::ComputeFeatures(const vector<BaseFloat> &wave, BaseFloat sample_freq, BaseFloat vtln_warp, vector<BaseFloat* > &output)
{
    int32 rows_out = NumFrames(wave.size());
    int32 cols_out = mfccOptions.num_ceps;

    output.resize(rows_out, NULL);

    vector<BaseFloat> window;  // windowed waveform.
    for (int32 frame = 0; frame < rows_out; ++frame)
    {
        ExtractWindow(wave, frame, vtln_warp, window, output[frame]);
    }
}

int32 MfccComputer::NumFrames(int64 num_samples)
{
    int64 frame_shift = frameOptions.WindowShift();
    int64 frame_length = frameOptions.WindowSize();

    if (num_samples < frame_length)
    {
        return 0;
    }
    else
    {
        return (1 + ((num_samples - frame_length) / frame_shift));
    }
}

void MfccComputer::ExtractWindow(const vector<BaseFloat> &wave, int32 f, BaseFloat vtln_warp, vector<BaseFloat> &window, BaseFloat* output)
{
    int32 frame_length = frameOptions.WindowSize();
    int32 frame_length_padded = frameOptions.PaddedWindowSize();

    if(window.size() != frame_length_padded)
    {
        window.resize(frame_length_padded);
    }

    memcpy(window.data(), wave.data()+f*frameOptions.WindowShift(), frame_length*sizeof(BaseFloat));
    memset(window.data()+frame_length, 0, (frame_length_padded-frame_length)*sizeof(BaseFloat));

    ProcessWindow(window, vtln_warp, output);
}

void MfccComputer::ProcessWindow(vector<BaseFloat> window, BaseFloat vtln_warp, BaseFloat* output)
{
    int32 frame_length = frameOptions.WindowSize();

    if (frameOptions.dither != 0.0)
    {
        Dither(window, frame_length, frameOptions.dither);
    }

    if (frameOptions.remove_dc_offset)
    {
        BaseFloat offset = -Sum(window) / frame_length;
        for(int i=0; i<frame_length; i++)
        {
            window[i] += offset;
        }
    }

    if (frameOptions.preemph_coeff != 0.0)
    {
        Preemphasize(window, frame_length, frameOptions.preemph_coeff);
    }

    for(int i=0; i<frame_length; i++)
    {
        window[i] *= frameOptions.window[i];
    }

    srfft->Compute(window.data(), true);

    ComputePowerSpectrum(window);

    const MelBanks &mel_banks = *(GetMelBanks(vtln_warp));
    mel_banks.Compute(window, mel_energies_);

    ApplyFloor(mel_energies_, std::numeric_limits<float>::epsilon());
    ApplyLog(mel_energies_);

    output = (BaseFloat*)malloc(mfccOptions.num_ceps * sizeof(BaseFloat));

    for(int32 i=0; i<mfccOptions.num_ceps; i++)
    {
        output[i] = 0.0f;
        for(int32 j=0; j<mel_energies_.size(); j++)
        {
            output[i] += mel_energies_[j] * dct_matrix.data[i*dct_matrix.cols + j];
        }
    }

    if (mfccOptions.cepstral_lifter != 0.0)
    {
        for(int32 i=0; i<mfccOptions.num_ceps; i++)
        {
            output[i] *= lifter_coeffs_[i];
        }
    }
}

const MelBanks *MfccComputer::GetMelBanks(BaseFloat vtln_warp) {
  MelBanks *this_mel_banks = NULL;
  std::map<BaseFloat, MelBanks*>::iterator iter = mel_banks_.find(vtln_warp);
  if (iter == mel_banks_.end()) {
    this_mel_banks = new MelBanks(mfccOptions.mel_opts,
                                  frameOptions,
                                  vtln_warp);
    mel_banks_[vtln_warp] = this_mel_banks;
  } else {
    this_mel_banks = iter->second;
  }
  return this_mel_banks;
}