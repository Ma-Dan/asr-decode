#ifndef SIMPLE_DECODER
#define SIMPLE_DECODER

#include "common.h"
#include "transition-model.h"
#include "am-diag-gmm.h"
#include "fstreader.h"

typedef int StateId;

typedef struct tagDecodeArc
{
    int ilabel;
    int olabel;
    BaseFloat weight1;
    BaseFloat weight2;
    int nextstate;
} DecodeArc, *P_DecodeArc;

typedef struct Token
{
    DecodeArc arc;
    Token *prev;
    int32 ref_count;
    double cost;
} *P_Token;

class SimpleDecoder
{
    public:
        SimpleDecoder(TransitionModel *transmodel, AmDiagGmm *amgmm, FstReader *fst, BaseFloat beam);
        bool Decode(P_Matrix feature, BaseFloat acoustic_scale);
        vector<int> GetBestPath();

        void InitDecoding();
        void AdvanceDecoding(P_Matrix feature, BaseFloat acoustic_scale);

    private:
        class TransitionModel *m_transmodel;
        class AmDiagGmm *m_amgmm;
        map<StateId, Token*> cur_toks;
        map<StateId, Token*> prev_toks;
        class FstReader *m_fst;
        BaseFloat m_beam;
        int32 num_frames_decoded;

        void ProcessEmitting(P_Matrix feature, BaseFloat acoustic_scale);
        void ProcessNonemitting();

        BaseFloat LogLikelihood(P_Matrix feature, int32 frame, int32 tid);

        static void ClearToks(map<StateId, Token*> &toks);
        static void PruneToks(BaseFloat beam, map<StateId, Token*> *toks);
};

#endif