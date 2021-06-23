#include "fstreader.h"

int main(int argc, char* argv[])
{
    const char* fstFileName = "/Users/dan/Documents/ASR/kaldi/egs/yesno/s5/exp/mono0a/graph_tgpr/HCLG.fst";
    //const char* fstFileName = "/Users/dan/Documents/ASR/vosk-api/c/model/HCLG.fst";
    FstReader fstReader;
    fstReader.Read(fstFileName);
    return 0;
}