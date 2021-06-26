#include "transition-model.h"

const char* transmodel = "<TransitionModel>";
const char* topology = "<Topology>";
const char* tuplesName = "<Tuples>";
const char* triplesName = "<Triples>";

void TransitionModel::Read(FILE* fp)
{
    char token[128];
    ReadToken(fp, token);
    if(strcmp(transmodel, token) != 0)
    {
        printf("Model file type error!\n");
        return;
    }

    ReadTopo(fp);

    //Read tuples
    ReadToken(fp, token);
    int32 size;
    ReadBasicType(fp, &size);
    tuples.resize(size);
    for (int32 i = 0; i < size; i++)
    {
        ReadBasicType(fp, &(tuples[i].phone));
        ReadBasicType(fp, &(tuples[i].hmm_state));
        ReadBasicType(fp, &(tuples[i].forward_pdf));
        if (0 == strcmp(token, tuplesName))
        {
            ReadBasicType(fp, &(tuples[i].self_loop_pdf));
        }
        else if (0 == strcmp(token, triplesName))
        {
            tuples[i].self_loop_pdf = tuples[i].forward_pdf;
        }
    }
    ReadToken(fp, token);
    //TODO: Check token is </Triples> or </Tuples>
    ComputeDerived();
    ReadToken(fp, token); //<LogProbs>
    ReadFloatVectors(fp, &log_probs);
    ReadToken(fp, token); //</LogProbs>
    ReadToken(fp, token); //</TransitionModel>
    ComputeDerivedOfProbs();
    //TODO: Check
}

void TransitionModel::ReadTopo(FILE *fp)
{
    char token[128];
    ReadToken(fp, token);
    if(strcmp(topology, token) != 0)
    {
        printf("Topology file type error!\n");
        return;
    }

    ReadIntegerVector(fp, &topo.phones);
    ReadIntegerVector(fp, &topo.phone2idx);

    //Read Tuples
    int32 size;
    ReadBasicType(fp, &size);
    bool is_hmm = true;
    topo.entries.resize(size);
    for (int32 i = 0; i < size; i++)
    {
        int32 thist_sz;
        ReadBasicType(fp, &thist_sz);
        topo.entries[i].resize(thist_sz);
        for (int32 j = 0 ; j < thist_sz; j++)
        {
            ReadBasicType(fp, &(topo.entries[i][j].forward_pdf_class));
            if(is_hmm)
            {
                topo.entries[i][j].self_loop_pdf_class = topo.entries[i][j].forward_pdf_class;
            }
            else
            {
                ReadBasicType(fp, &(topo.entries[i][j].self_loop_pdf_class));
            }
            int32 thiss_sz;
            ReadBasicType(fp, &thiss_sz);
            topo.entries[i][j].transitions.resize(thiss_sz);
            for (int32 k = 0; k < thiss_sz; k++)
            {
                ReadBasicType(fp, &(topo.entries[i][j].transitions[k].first));
                ReadBasicType(fp, &(topo.entries[i][j].transitions[k].second));
            }
        }
    }
    ReadToken(fp, token);
    //TODO: Add check
}

void TransitionModel::ComputeDerived()
{
    state2id.resize(tuples.size()+2);  // indexed by transition-state, which
    // is one based, but also an entry for one past end of list.

    int32 cur_transition_id = 1;
    num_pdfs = 0;
    for (int32 tstate = 1;
        tstate <= static_cast<int32>(tuples.size()+1);  // not a typo.
        tstate++)
    {
        state2id[tstate] = cur_transition_id;
        if (static_cast<size_t>(tstate) <= tuples.size())
        {
          int32 phone = tuples[tstate-1].phone,
              hmm_state = tuples[tstate-1].hmm_state,
              forward_pdf = tuples[tstate-1].forward_pdf,
              self_loop_pdf = tuples[tstate-1].self_loop_pdf;
          num_pdfs = max(num_pdfs, 1 + forward_pdf);
          num_pdfs = max(num_pdfs, 1 + self_loop_pdf);
          const HmmState &state = TopologyForPhone(phone)[hmm_state];
          int32 my_num_ids = static_cast<int32>(state.transitions.size());
          cur_transition_id += my_num_ids;  // # trans out of this state.
        }
    }

    id2state.resize(cur_transition_id);   // cur_transition_id is #transition-ids+1.
    id2pdf_id.resize(cur_transition_id);
    for (int32 tstate = 1; tstate <= static_cast<int32>(tuples.size()); tstate++)
    {
        for (int32 tid = state2id[tstate]; tid < state2id[tstate+1]; tid++)
        {
            id2state[tid] = tstate;
            if (IsSelfLoop(tid))
            {
                id2pdf_id[tid] = tuples[tstate-1].self_loop_pdf;
            }
            else
            {
                id2pdf_id[tid] = tuples[tstate-1].forward_pdf;
            }
        }
    }

    // The following statements put copies a large number in the region of memory
    // past the end of the id2pdf_id_ array, while leaving the array as it was
    // before.  The goal of this is to speed up decoding by disabling a check
    // inside TransitionIdToPdf() that the transition-id was within the correct
    // range.
    int32 num_big_numbers = min<int32>(2000, cur_transition_id);
    id2pdf_id.resize(cur_transition_id + num_big_numbers,
                      std::numeric_limits<int32>::max());
    id2pdf_id.resize(cur_transition_id);
}

bool TransitionModel::IsSelfLoop(int32 trans_id) const {
    int32 trans_state = id2state[trans_id];
    int32 trans_index = trans_id - state2id[trans_state];
    const Tuple &tuple = tuples[trans_state-1];
    int32 phone = tuple.phone, hmm_state = tuple.hmm_state;
    const TopologyEntry &entry = TopologyForPhone(phone);
    return (static_cast<size_t>(trans_index) < entry[hmm_state].transitions.size()
            && entry[hmm_state].transitions[trans_index].first == hmm_state);
}

const TopologyEntry& TransitionModel::TopologyForPhone(int32 phone) const
{
    // Will throw if phone not covered.
    if (static_cast<size_t>(phone) >= topo.phone2idx.size() || topo.phone2idx[phone] == -1) {
        printf("TopologyForPhone(), phone %d not covered.\n", phone);
    }
    return topo.entries[topo.phone2idx[phone]];
}

void TransitionModel::ComputeDerivedOfProbs()
{
    non_self_loop_log_probs.resize(NumTransitionStates()+1);  // this array indexed
    //  by transition-state with nothing in zeroth element.
    for (int32 tstate = 1; tstate <= NumTransitionStates(); tstate++)
    {
        int32 tid = SelfLoopOf(tstate);
        if (tid == 0)
        {   // no self-loop
            non_self_loop_log_probs[tstate] = 0.0;  // log(1.0)
        }
        else
        {
            float self_loop_prob = expf(GetTransitionLogProb(tid)),
                  non_self_loop_prob = 1.0 - self_loop_prob;
            if (non_self_loop_prob <= 0.0)
            {
                printf("ComputeDerivedOfProbs(): non-self-loop prob is %f\n", non_self_loop_prob);
                non_self_loop_prob = 1.0e-10;  // just so we can continue...
            }
            non_self_loop_log_probs[tstate] = logf(non_self_loop_prob);  // will be negative.
        }
  }
}

int32 TransitionModel::SelfLoopOf(int32 trans_state) const
{   // returns the self-loop transition-id
    const Tuple &tuple = tuples[trans_state-1];
    // or zero if does not exist.
    int32 phone = tuple.phone, hmm_state = tuple.hmm_state;
    const TopologyEntry &entry = TopologyForPhone(phone);

    for (int32 trans_index = 0;
        trans_index < static_cast<int32>(entry[hmm_state].transitions.size());
        trans_index++)
    {
        if (entry[hmm_state].transitions[trans_index].first == hmm_state)
        {
            return PairToTransitionId(trans_state, trans_index);
        }
    }

    return 0;  // invalid transition id.
}

float TransitionModel::GetTransitionLogProb(int32 trans_id) const
{
  return log_probs[trans_id];
}

int32 TransitionModel::PairToTransitionId(int32 trans_state, int32 trans_index) const
{
    return state2id[trans_state] + trans_index;
}

int32 TransitionModel::NumTransitionStates()
{
    return tuples.size();
}