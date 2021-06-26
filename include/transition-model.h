#ifndef TRANSITION_MODEL
#define TRANSITION_MODEL

#include "common.h"

typedef struct tagHmmState
{
    int32 forward_pdf_class;
    int32 self_loop_pdf_class;
    vector<pair<int32, float> > transitions;
} HmmState, *P_HmmState;

typedef vector<HmmState> TopologyEntry;

typedef struct tagHmmTopology
{
    vector<int32> phones;
    vector<int32> phone2idx;
    vector<TopologyEntry> entries;
} HmmTopology, *P_HmmTopology;

typedef struct tagTuple
{
    int32 phone;
    int32 hmm_state;
    int32 forward_pdf;
    int32 self_loop_pdf;
} Tuple, *P_Tuple;

class TransitionModel
{
    public:
        void Read(FILE *fp);
        int32 TransitionIdToPdf(int32 trans_id) const;

    private:
        HmmTopology topo;
        vector<Tuple> tuples;
        vector<int32> state2id;
        vector<int32> id2state;
        vector<int32> id2pdf_id;
        vector<float> log_probs;
        vector<float> non_self_loop_log_probs;
        int32 num_pdfs;

        void ReadTopo(FILE *fp);

        void ComputeDerived();
        bool IsSelfLoop(int32 trans_id) const;
        const TopologyEntry& TopologyForPhone(int32 phone) const;
        void ComputeDerivedOfProbs();
        int32 SelfLoopOf(int32 trans_state) const;
        float GetTransitionLogProb(int32 trans_id) const;
        int32 PairToTransitionId(int32 trans_state, int32 trans_index) const;
        int32 NumTransitionStates();
};

#endif