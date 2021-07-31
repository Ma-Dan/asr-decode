#ifndef COMPRESSED_MATRIX
#define COMPRESSED_MATRIX

#include "common.h"

struct GlobalHeader {
    int32 format;     // Represents the enum DataFormat.
    float min_value;  // min_value and range represent the ranges of the integer
                      // data in the kTwoByte and kOneByte formats, and the
                      // range of the PerColHeader uint16's in the
                      // kOneByteWithColheaders format.
    float range;
    int32 num_rows;
    int32 num_cols;
};

struct PerColHeader {
    uint16 percentile_0;
    uint16 percentile_25;
    uint16 percentile_75;
    uint16 percentile_100;
};

class CompressedMatrix {
    public:
        CompressedMatrix(): data_(NULL) { }
        ~CompressedMatrix() { Clear(); }

        void Clear();

        void CopyFromMat(const P_Matrix mat);
        void ComputeGlobalHeader(const P_Matrix mat, GlobalHeader *header);
        static int32 DataSize(const GlobalHeader &header);
        static void* AllocateData(int32 num_bytes);
        void CopyToMat(P_Matrix mat) const;

    private:
        void GetMinMax(const P_Matrix mat, BaseFloat* pMin, BaseFloat* pMax);
        static void CompressColumn(const GlobalHeader &global_header,
                             const BaseFloat *data, int32 stride,
                             int32 num_rows, PerColHeader *header,
                             uint8 *byte_data);
        static void ComputeColHeader(const GlobalHeader &global_header,
                               const BaseFloat *data, int32 stride,
                               int32 num_rows, PerColHeader *header);

        static inline uint16 FloatToUint16(const GlobalHeader &global_header,
                                     float value);

        static inline float Uint16ToFloat(const GlobalHeader &global_header,
                                    uint16 value);

        // this is used only in the kOneByteWithColHeaders compression format.
        static inline uint8 FloatToChar(float p0, float p25,
                                          float p75, float p100,
                                          float value);

        // this is used only in the kOneByteWithColHeaders compression format.
        static inline float CharToFloat(float p0, float p25,
                                  float p75, float p100,
                                  uint8 value);

        void *data_;
};

#endif