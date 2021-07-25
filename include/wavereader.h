#ifndef WAVE_READER
#define WAVE_READER

#include "common.h"

typedef struct tagWaveHeader
{
    uint8  chunk_id[4];      //'RIFF'
    uint32 chunk_size;
    uint8  format[4];        //'WAVE'
    uint8  subchunk1_id[4];  //'FMT'
    uint32 subchunk1_size;   //PCM = 16
    uint16 audio_format;     //PCM = 1
    uint16 channels;
    uint32 sample_rate;
    uint32 byte_rate;
    uint16 block_align;      //NumChannels * BitsPerSample / 8
    uint16 bit_per_sample;
    uint8  subchunk2_id[4];  //'DATA'
    uint32 subchunk2_size;
} WaveHeader, *P_WaveHeader;

typedef struct tagWaveFile
{
    WaveHeader header;
    void       *data;
} WaveFile, *P_WaveFile;

class WaveReader
{
    public:
        WaveReader();
        ~WaveReader();
        void ReadWaveFile(const char* fileName);

    private:
        WaveFile m_wavefile;
};

#endif