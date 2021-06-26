#include "transition-model.h"
#include "am-diag-gmm.h"
#include "fstreader.h"
#include "simple-decoder.h"

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        printf("arg error\n");
        return -1;
    }

    char* mdlFileName = argv[1];
    char* fstFileName = argv[2];

    float beam = 16.0;

    // Read Transition model and GMM AM model
    FILE *fpMdl = fopen(mdlFileName, "rb");

    bool binary = false;
    char hdr[2];
    fread(hdr, 2, 1, fpMdl);
    if(hdr[1] == 'B')
    {
        binary = true;
    }
    TransitionModel trans_model;
    trans_model.Read(fpMdl);

    AmDiagGmm am_gmm;
    am_gmm.Read(fpMdl);

    fclose(fpMdl);

    // Read HCLG fst
    FstReader fstReader;
    fstReader.Read(fstFileName);

    // Decode feature
    SimpleDecoder decoder(&fstReader, beam);
    decoder.Decode();

    return 0;
}