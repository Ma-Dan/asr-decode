#include "wavereader.h"

WaveReader::WaveReader()
{
    memset(&m_wavefile, 0, sizeof(WaveFile));
}

WaveReader::~WaveReader()
{
    m_wavefile.data.clear();
    m_waveData.clear();
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
    m_wavefile.data.resize(dataSize/2);
    fread(m_wavefile.data.data(), dataSize, 1, fp);

    m_waveData.clear();
    for(int i=0; i<dataSize/2; i++)
    {
        m_waveData.push_back(static_cast<BaseFloat>(m_wavefile.data[i]));
    }

    fclose(fp);
}