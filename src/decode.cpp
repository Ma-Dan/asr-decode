#include "wavereader.h"
#include "transition-model.h"
#include "am-diag-gmm.h"
#include "fstreader.h"
#include "feature-mfcc.h"
#include "compressed-matrix.h"
#include "compute-cmvn-stats.h"
#include "add-deltas.h"
#include "simple-decoder.h"

int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        printf("arg error\n");
        return -1;
    }

    char* mdlFileName = argv[1];
    char* fstFileName = argv[2];
    char* waveFileName = argv[3];

    BaseFloat vtln_warp = 1.0;

    BaseFloat acoustic_scale = 0.083333;
    BaseFloat beam = 16.0;

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

    // Read wave file
    WaveReader waveReader;
    waveReader.ReadWaveFile(waveFileName);

    // Compute MFCC
    MfccComputer mfccComputer;
    Matrix feats;
    mfccComputer.ComputeFeatures(waveReader.m_waveData, waveReader.m_wavefile.header.sample_rate, vtln_warp, &feats);

    // Compress matrix
    CompressedMatrix compressedMatrix;
    compressedMatrix.CopyFromMat(&feats);
    compressedMatrix.CopyToMat(&feats);

    // Compute CMVN stats and apply
    MatrixDouble cmvn_stats;
    InitCmvnStats(feats.cols, &cmvn_stats);
    AccCmvnStats(&feats, &cmvn_stats);
    ApplyCmvn(&cmvn_stats, false, &feats);

    // Add deltas
    DeltaFeaturesOptions opts;
    Matrix feature;
    ComputeDeltas(opts, &feats, &feature);

    // Decode feature
    SimpleDecoder decoder(&trans_model, &am_gmm, &fstReader, beam);
    decoder.Decode(&feature, acoustic_scale);

    vector<int> result = decoder.GetBestPath();

    printf("Decoded result: ");
    for(int i=0; i<result.size(); i++)
    {
        printf("%d ", result[i]);
    }
    printf("\n");

    return 0;
}