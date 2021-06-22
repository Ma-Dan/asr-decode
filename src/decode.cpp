#include "fstreader.h"

int main(int argc, char* argv[])
{
    const char* fstFileName = "/Users/dan/Documents/ASR/kaldi/egs/yesno/s5/exp/mono0a/graph_tgpr/HCLG.fst";
    FstHeader hdr;
    hdr.Read(fstFileName);
    return 0;
}