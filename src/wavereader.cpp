#include "wavereader.h"

WaveReader::WaveReader()
{
    memset(&m_wavefile, 0, sizeof(WaveFile));
}

WaveReader::~WaveReader()
{
    SAFE_FREE(m_wavefile.data);
}

void WaveReader::ReadWaveFile(const char* fileName)
{
    FILE *fp = fopen(fileName, "rb");

    if(!fp)
    {
        printf("Open wave file %s error\n", fileName);
    }

    //读取文件头
    fread(&m_wavefile, sizeof(WaveHeader), 1, fp);

    //读取数据
    int dataSize = m_wavefile.header.subchunk2_size;
    m_wavefile.data = malloc(dataSize);
    fread(&m_wavefile.data, dataSize, 1, fp);

    fclose(fp);
}