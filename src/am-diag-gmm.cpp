#include "am-diag-gmm.h"

void AmDiagGmm::Read(FILE *fp)
{
    int32 num_pdfs, dim;
    char token[128];

    ReadToken(fp, token); //<DIMENSION>
    ReadBasicType(fp, &dim);
    ReadToken(fp, token); //<NUMPDFS>
    ReadBasicType(fp, &num_pdfs);

    densities.reserve(num_pdfs);
    for (int32 i = 0; i < num_pdfs; i++)
    {
        densities.push_back(ReadDiagGmm(fp));
    }
}

AmDiagGmm::~AmDiagGmm()
{
}

DiagGmm* AmDiagGmm::ReadDiagGmm(FILE *fp)
{
    DiagGmm *diag_gmm = new DiagGmm();
    char token[128];

    ReadToken(fp, token); //<DiagGMMBegin> or <DiagGMM>

    ReadToken(fp, token);
    if(0 == strcmp(token, "<GCONSTS>"))
    {
        ReadFloatVectors(fp, &diag_gmm->gconsts);
    }

    ReadToken(fp, token);
    if(0 == strcmp(token, "<WEIGHTS>"))
    {
        ReadFloatVectors(fp, &diag_gmm->weights);
    }

    ReadToken(fp, token); //<MEANS_INVVARS>
    ReadFloatMatrix(fp, &diag_gmm->means_invvars);

    ReadToken(fp, token); //<INV_VARS>
    ReadFloatMatrix(fp, &diag_gmm->inv_vars);

    ReadToken(fp, token); //</DiagGMM>

    ComputeGconsts(diag_gmm);

    return diag_gmm;
}

int32 AmDiagGmm::ComputeGconsts(DiagGmm* diaggmm) {
    int32 num_mix = NumGauss(diaggmm);
    int32 dim = Dim(diaggmm);
    float offset = -0.5 * M_LOG_2PI * dim;  // constant term in gconst.
    int32 num_bad = 0;

    // Resize if Gaussians have been removed during Update()
    if (num_mix != static_cast<int32>(diaggmm->gconsts.size()))
    {
        diaggmm->gconsts.resize(num_mix);
    }

    for (int32 mix = 0; mix < num_mix; mix++)
    {
        float gc = logf(diaggmm->weights[mix]) + offset;  // May be -inf if weights == 0
        for (int32 d = 0; d < dim; d++)
        {
            gc += 0.5 * logf(ReadMatrix(&diaggmm->inv_vars, mix, d)) - 0.5 * ReadMatrix(&diaggmm->means_invvars, mix, d)
                  * ReadMatrix(&diaggmm->means_invvars, mix, d) / ReadMatrix(&diaggmm->inv_vars, mix, d);
        }
        // Change sign for logdet because var is inverted. Also, note that
        // mean_invvars(mix, d)*mean_invvars(mix, d)/inv_vars(mix, d) is the
        // mean-squared times inverse variance, since mean_invvars(mix, d) contains
        // the mean times inverse variance.
        // So gc is the likelihood at zero feature value.

        if (isnan(gc))
        {  // negative infinity is OK but NaN is not acceptable
           printf("At component %d  not a number in gconst computation", mix);
        }
        if (isinf(gc))
        {
            num_bad++;
            // If positive infinity, make it negative infinity.
            // Want to make sure the answer becomes -inf in the end, not NaN.
            if (gc > 0)
            {
                gc = -gc;
            }
        }
        diaggmm->gconsts[mix] = gc;
    }

    diaggmm->valid_gconsts = true;
    return num_bad;
}

int32 AmDiagGmm::NumGauss(DiagGmm* diaggmm) const
{
    return diaggmm->weights.size();
}

int32 AmDiagGmm::Dim(DiagGmm* diaggmm) const
{
    return diaggmm->means_invvars.cols;
}