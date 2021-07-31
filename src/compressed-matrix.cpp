#include "compressed-matrix.h"

void CompressedMatrix::Clear() {
  if (data_ != NULL) {
    delete [] static_cast<float*>(data_);
    data_ = NULL;
  }
}

int32 CompressedMatrix::DataSize(const GlobalHeader &header) {
  return sizeof(GlobalHeader) +
        header.num_cols * (sizeof(PerColHeader) + header.num_rows);
}

void* CompressedMatrix::AllocateData(int32 num_bytes) {
  return reinterpret_cast<void*>(new float[(num_bytes/3) + 4]);
}

void CompressedMatrix::CopyFromMat(const P_Matrix mat)
{
    Clear();

    GlobalHeader global_header;
    ComputeGlobalHeader(mat, &global_header);

    int32 data_size = DataSize(global_header);

    data_ = AllocateData(data_size);

    *(reinterpret_cast<GlobalHeader*>(data_)) = global_header;

    PerColHeader *header_data =
        reinterpret_cast<PerColHeader*>(static_cast<char*>(data_) +
                                        sizeof(GlobalHeader));
    uint8 *byte_data =
        reinterpret_cast<uint8*>(header_data + global_header.num_cols);

    const BaseFloat *matrix_data = mat->data.data();

    for (int32 col = 0; col < global_header.num_cols; col++) {
      CompressColumn(global_header,
                     matrix_data + col, mat->cols,
                     global_header.num_rows,
                     header_data, byte_data);
      header_data++;
      byte_data += global_header.num_rows;
    }
}

void CompressedMatrix::CopyToMat(P_Matrix mat) const {
  GlobalHeader *h = reinterpret_cast<GlobalHeader*>(data_);
  int32 num_cols = h->num_cols, num_rows = h->num_rows;
  if (1) {
    PerColHeader *per_col_header = reinterpret_cast<PerColHeader*>(h+1);
    uint8 *byte_data = reinterpret_cast<uint8*>(per_col_header +
                                                h->num_cols);
    for (int32 i = 0; i < num_cols; i++, per_col_header++) {
      float p0 = Uint16ToFloat(*h, per_col_header->percentile_0),
          p25 = Uint16ToFloat(*h, per_col_header->percentile_25),
          p75 = Uint16ToFloat(*h, per_col_header->percentile_75),
          p100 = Uint16ToFloat(*h, per_col_header->percentile_100);
      for (int32 j = 0; j < num_rows; j++, byte_data++) {
        float f = CharToFloat(p0, p25, p75, p100, *byte_data);
        mat->data[j*num_cols+i] = f;
      }
    }
  }
}

void CompressedMatrix::ComputeGlobalHeader(const P_Matrix mat, GlobalHeader *header)
{
    header->num_rows = mat->rows;
    header->num_cols = mat->cols;

    BaseFloat min_value, max_value;
    GetMinMax(mat, &min_value, &max_value);

    header->min_value = min_value;
    header->range = max_value - min_value;
}

void CompressedMatrix::GetMinMax(const P_Matrix mat, BaseFloat* pMin, BaseFloat* pMax)
{
    int32 total = mat->rows * mat->cols;

    *pMin = mat->data[0];
    *pMax = mat->data[0];

    for(int32 i=1; i<total; i++)
    {
        if(*pMin > mat->data[i])
        {
            *pMin = mat->data[i];
        }

        if(*pMax < mat->data[i])
        {
            *pMax = mat->data[i];
        }
    }
}

void CompressedMatrix::CompressColumn(
    const GlobalHeader &global_header,
    const BaseFloat *data, int32 stride,
    int32 num_rows, PerColHeader *header,
    uint8 *byte_data) {
  ComputeColHeader(global_header, data, stride,
                   num_rows, header);

  float p0 = Uint16ToFloat(global_header, header->percentile_0),
      p25 = Uint16ToFloat(global_header, header->percentile_25),
      p75 = Uint16ToFloat(global_header, header->percentile_75),
      p100 = Uint16ToFloat(global_header, header->percentile_100);

  for (int32 i = 0; i < num_rows; i++) {
    BaseFloat this_data = data[i * stride];
    byte_data[i] = FloatToChar(p0, p25, p75, p100, this_data);
  }
}

void CompressedMatrix::ComputeColHeader(
    const GlobalHeader &global_header,
    const BaseFloat *data, int32 stride,
    int32 num_rows, PerColHeader *header) {
  std::vector<BaseFloat> sdata(num_rows); // the sorted data.
  for (size_t i = 0, size = sdata.size(); i < size; i++)
    sdata[i] = data[i*stride];

  if (num_rows >= 5) {
    int quarter_nr = num_rows/4;
    // std::sort(sdata.begin(), sdata.end());
    // The elements at positions 0, quarter_nr,
    // 3*quarter_nr, and num_rows-1 need to be in sorted order.
    std::nth_element(sdata.begin(), sdata.begin() + quarter_nr, sdata.end());
    // Now, sdata.begin() + quarter_nr contains the element that would appear
    // in sorted order, in that position.
    std::nth_element(sdata.begin(), sdata.begin(), sdata.begin() + quarter_nr);
    // Now, sdata.begin() and sdata.begin() + quarter_nr contain the elements
    // that would appear at those positions in sorted order.
    std::nth_element(sdata.begin() + quarter_nr + 1,
                     sdata.begin() + (3*quarter_nr), sdata.end());
    // Now, sdata.begin(), sdata.begin() + quarter_nr, and sdata.begin() +
    // 3*quarter_nr, contain the elements that would appear at those positions
    // in sorted order.
    std::nth_element(sdata.begin() + (3*quarter_nr) + 1, sdata.end() - 1,
                     sdata.end());
    // Now, sdata.begin(), sdata.begin() + quarter_nr, and sdata.begin() +
    // 3*quarter_nr, and sdata.end() - 1, contain the elements that would appear
    // at those positions in sorted order.

    header->percentile_0 =
        std::min<uint16>(FloatToUint16(global_header, sdata[0]), 65532);
    header->percentile_25 =
        std::min<uint16>(
            std::max<uint16>(
                FloatToUint16(global_header, sdata[quarter_nr]),
                header->percentile_0 + static_cast<uint16>(1)), 65533);
    header->percentile_75 =
        std::min<uint16>(
            std::max<uint16>(
                FloatToUint16(global_header, sdata[3*quarter_nr]),
                header->percentile_25 + static_cast<uint16>(1)), 65534);
    header->percentile_100 = std::max<uint16>(
        FloatToUint16(global_header, sdata[num_rows-1]),
        header->percentile_75 + static_cast<uint16>(1));

  }
}

inline uint16 CompressedMatrix::FloatToUint16(
    const GlobalHeader &global_header,
    float value) {
  float f = (value - global_header.min_value) /
      global_header.range;
  if (f > 1.0) f = 1.0;  // Note: this should not happen.
  if (f < 0.0) f = 0.0;  // Note: this should not happen.
  return static_cast<int>(f * 65535 + 0.499);  // + 0.499 is to
  // round to closest int; avoids bias.
}

inline uint8 CompressedMatrix::FloatToChar(
    float p0, float p25, float p75, float p100,
    float value) {
  int ans;
  if (value < p25) {  // range [ p0, p25 ) covered by
    // characters 0 .. 64.  We round to the closest int.
    float f = (value - p0) / (p25 - p0);
    ans = static_cast<int>(f * 64 + 0.5);
    // Note: the checks on the next two lines
    // are necessary in pathological cases when all the elements in a row
    // are the same and the percentile_* values are separated by one.
    if (ans < 0) ans = 0;
    if (ans > 64) ans = 64;
  } else if (value < p75) {  // range [ p25, p75 )covered
    // by characters 64 .. 192.  We round to the closest int.
    float f = (value - p25) / (p75 - p25);
    ans = 64 + static_cast<int>(f * 128 + 0.5);
    if (ans < 64) ans = 64;
    if (ans > 192) ans = 192;
  } else {  // range [ p75, p100 ] covered by
    // characters 192 .. 255.  Note: this last range
    // has fewer characters than the left range, because
    // we go up to 255, not 256.
    float f = (value - p75) / (p100 - p75);
    ans = 192 + static_cast<int>(f * 63 + 0.5);
    if (ans < 192) ans = 192;
    if (ans > 255) ans = 255;
  }
  return static_cast<uint8>(ans);
}

inline float CompressedMatrix::Uint16ToFloat(
    const GlobalHeader &global_header,
    uint16 value) {
  // the constant 1.52590218966964e-05 is 1/65535.
  return global_header.min_value
      + global_header.range * 1.52590218966964e-05F * value;
}

inline float CompressedMatrix::CharToFloat(
    float p0, float p25, float p75, float p100,
    uint8 value) {
  if (value <= 64) {
    return p0 + (p25 - p0) * value * (1/64.0);
  } else if (value <= 192) {
    return p25 + (p75 - p25) * (value - 64) * (1/128.0);
  } else {
    return p75 + (p100 - p75) * (value - 192) * (1/63.0);
  }
}
