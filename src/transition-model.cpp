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

void TransitionModel::ReadIntegerVector(FILE *fp, vector<int32> *v)
{
    uint8 size = 0;
    fread(&size, sizeof(size), 1, fp);

    if(size != sizeof(int32))
    {
        printf("vector size error!\n");
        return;
    }

    uint32 vsize = 0;
    fread(&vsize, sizeof(vsize), 1, fp);

    int32 value;
    for(int i=0; i<vsize; i++)
    {
        fread(&value, sizeof(value), 1, fp);
        v->push_back(value);
    }
}

void TransitionModel::ReadBasicType(FILE *fp, int32 *t)
{
    uint8 size = 0;
    fread(&size, sizeof(size), 1, fp);

    if(size != sizeof(int32))
    {
        printf("int32 size error!\n");
        return;
    }

    fread(t, sizeof(*t), 1, fp);
}

void TransitionModel::ReadBasicType(FILE *fp, float *t)
{
    uint8 size = 0;
    fread(&size, sizeof(size), 1, fp);

    if(size != sizeof(float))
    {
        printf("float size error!\n");
        return;
    }

    fread(t, sizeof(*t), 1, fp);
}

void TransitionModel::ReadToken(FILE *fp, char* s)
{
    int index = 0;
    char c = '\0';
    while(c != ' ')
    {
        fread(&c, 1, 1, fp);
        s[index] = c;
        index++;
    }

    s[index-1] = '\0';
}

/*void TransitionModel::ComputeDerived()
{
    state2id.resize(tuples.size()+2);  // indexed by transition-state, which
    // is one based, but also an entry for one past end of list.

    int32 cur_transition_id = 1;
    num_pdfs = 0;
  for (int32 tstate = 1;
      tstate <= static_cast<int32>(tuples_.size()+1);  // not a typo.
      tstate++) {
    state2id_[tstate] = cur_transition_id;
    if (static_cast<size_t>(tstate) <= tuples_.size()) {
      int32 phone = tuples_[tstate-1].phone,
          hmm_state = tuples_[tstate-1].hmm_state,
          forward_pdf = tuples_[tstate-1].forward_pdf,
          self_loop_pdf = tuples_[tstate-1].self_loop_pdf;
      num_pdfs_ = std::max(num_pdfs_, 1 + forward_pdf);
      num_pdfs_ = std::max(num_pdfs_, 1 + self_loop_pdf);
      const HmmTopology::HmmState &state = topo_.TopologyForPhone(phone)[hmm_state];
      int32 my_num_ids = static_cast<int32>(state.transitions.size());
      cur_transition_id += my_num_ids;  // # trans out of this state.
    }
  }

  id2state_.resize(cur_transition_id);   // cur_transition_id is #transition-ids+1.
  id2pdf_id_.resize(cur_transition_id);
  for (int32 tstate = 1; tstate <= static_cast<int32>(tuples_.size()); tstate++) {
    for (int32 tid = state2id_[tstate]; tid < state2id_[tstate+1]; tid++) {
      id2state_[tid] = tstate;
      if (IsSelfLoop(tid))
        id2pdf_id_[tid] = tuples_[tstate-1].self_loop_pdf;
      else
        id2pdf_id_[tid] = tuples_[tstate-1].forward_pdf;
    }
  }

  // The following statements put copies a large number in the region of memory
  // past the end of the id2pdf_id_ array, while leaving the array as it was
  // before.  The goal of this is to speed up decoding by disabling a check
  // inside TransitionIdToPdf() that the transition-id was within the correct
  // range.
  int32 num_big_numbers = std::min<int32>(2000, cur_transition_id);
  id2pdf_id_.resize(cur_transition_id + num_big_numbers,
                    std::numeric_limits<int32>::max());
  id2pdf_id_.resize(cur_transition_id);
}

bool TransitionModel::IsSelfLoop(int32 trans_id) const {
  KALDI_ASSERT(static_cast<size_t>(trans_id) < id2state_.size());
  int32 trans_state = id2state_[trans_id];
  int32 trans_index = trans_id - state2id_[trans_state];
  const Tuple &tuple = tuples_[trans_state-1];
  int32 phone = tuple.phone, hmm_state = tuple.hmm_state;
  const HmmTopology::TopologyEntry &entry = topo_.TopologyForPhone(phone);
  KALDI_ASSERT(static_cast<size_t>(hmm_state) < entry.size());
  return (static_cast<size_t>(trans_index) < entry[hmm_state].transitions.size()
          && entry[hmm_state].transitions[trans_index].first == hmm_state);
}*/