#include "transition-model.h"
#include "am-diag-gmm.h"
#include "fstreader.h"
#include "simple-decoder.h"

void ReadFeature(const char* fileName, P_Matrix feature)
{
    FILE *fp = fopen(fileName, "rb");

    //Hard-coded to read one feature
    feature->rows = 668;
    feature->cols = 39;
    feature->stride = 40;

    feature->data.resize(feature->rows * feature->cols);
    fseek(fp, 31, SEEK_SET);

    fread(feature->data.data(), sizeof(float), feature->rows * feature->cols, fp);

    fclose(fp);
}

int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        printf("arg error\n");
        return -1;
    }

    char* mdlFileName = argv[1];
    char* fstFileName = argv[2];
    char* featureFileName = argv[3];

    float acoustic_scale = 0.083333;
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

    // Read feature
    Matrix feature;
    ReadFeature(featureFileName, &feature);

    // Decode feature
    SimpleDecoder decoder(&trans_model, &am_gmm, &fstReader, beam);
    decoder.Decode(&feature, acoustic_scale);

    return 0;
}