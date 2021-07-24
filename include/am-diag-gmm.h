#ifndef AM_DIAG_GMM
#define AM_DIAG_GMM

#include "common.h"

class DiagGmm
{
    public:
        vector<BaseFloat> gconsts;
        bool valid_gconsts;
        vector<BaseFloat> weights;
        Matrix inv_vars;
        Matrix means_invvars;
};

class AmDiagGmm {
    public:
        void Read(FILE *fp);
        ~AmDiagGmm();
        DiagGmm& GetPdf(int32 pdf_index) const;

    private:
        vector<DiagGmm*> densities;
        DiagGmm* ReadDiagGmm(FILE *fp);
        int32 ComputeGconsts(DiagGmm* diaggmm);
        int32 NumGauss(DiagGmm* diaggmm) const;
        int32 Dim(DiagGmm* diaggmm) const;
};

#endif