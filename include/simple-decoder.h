#ifndef SIMPLE_DECODER
#define SIMPLE_DECODER

#include "common.h"
#include "fstreader.h"

typedef int StateId;

typedef struct tagDecodeArc
{
    int ilabel;
    int olabel;
    float weight1;
    float weight2;
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
        SimpleDecoder(FstReader *fst, float beam);
        bool Decode();

        void InitDecoding();
        void AdvanceDecoding();

    private:
        map<StateId, Token*> cur_toks;
        map<StateId, Token*> prev_toks;
        class FstReader *m_fst;
        float m_beam;
        int32 num_frames_decoded;

        void ProcessEmitting();
        void ProcessNonemitting();

        static void ClearToks(map<StateId, Token*> &toks);
        static void PruneToks(float beam, map<StateId, Token*> *toks);
};

#endif