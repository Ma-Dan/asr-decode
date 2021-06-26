#include "simple-decoder.h"

static void TokenDelete(Token *tok)
{
    while(--tok->ref_count == 0)
    {
        Token *prev = tok->prev;
        SAFE_FREE(tok);
        if(prev == NULL)
        {
            return;
        }
        else tok = prev;
    }
}

static float GetDecodeArcWeight(P_DecodeArc arc)
{
    return arc->weight1 + arc->weight2;
}

static P_Token newToken(P_DecodeArc arc, float acoustic_cost, Token *prev)
{
    P_Token token = (P_Token)malloc(sizeof(Token));

    token->arc.ilabel = arc->ilabel;
    token->arc.olabel = arc->olabel;
    token->arc.weight1 = GetDecodeArcWeight(arc);
    token->arc.weight2 = acoustic_cost;
    token->arc.nextstate = arc->nextstate;

    token->prev = prev;
    token->ref_count = 1;

    if(prev)
    {
        prev->ref_count++;
        token->cost = prev->cost + (GetDecodeArcWeight(arc) + acoustic_cost);
    }
    else
    {
        token->cost = GetDecodeArcWeight(arc) + acoustic_cost;
    }

    return token;
}

static P_Token newToken(P_Arc arc, float acoustic_cost, Token *prev)
{
    P_Token token = (P_Token)malloc(sizeof(Token));

    token->arc.ilabel = arc->ilabel;
    token->arc.olabel = arc->olabel;
    token->arc.weight1 = arc->weight;
    token->arc.weight2 = acoustic_cost;
    token->arc.nextstate = arc->nextstate;

    token->prev = prev;
    token->ref_count = 1;

    if(prev)
    {
        prev->ref_count++;
        token->cost = prev->cost + (arc->weight + acoustic_cost);
    }
    else
    {
        token->cost = arc->weight + acoustic_cost;
    }

    return token;
}

SimpleDecoder::SimpleDecoder(TransitionModel *transmodel, AmDiagGmm *amgmm, FstReader *fst, float beam)
{
    m_transmodel = transmodel;
    m_amgmm = amgmm;
    m_fst = fst;
    m_beam = beam;
}

void SimpleDecoder::InitDecoding()
{
    // clean up from last time:
    ClearToks(cur_toks);
    ClearToks(prev_toks);
    // initialize decoding:
    StateId start_state = m_fst->Start();

    DecodeArc dummy_arc;
    dummy_arc.ilabel = 0;
    dummy_arc.olabel = 0;
    dummy_arc.weight1 = 0;
    dummy_arc.weight2 = 0;
    dummy_arc.nextstate = start_state;

    cur_toks[start_state] = newToken(&dummy_arc, 0.0, NULL);

    num_frames_decoded = 0;
    ProcessNonemitting();
}

bool SimpleDecoder::Decode(P_Matrix feature, float acoustic_scale)
{
    InitDecoding();
    AdvanceDecoding(feature, acoustic_scale);
    return (!cur_toks.empty());
}

vector<int> SimpleDecoder::GetBestPath()
{
    Token* best_token;
    float best_cost = std::numeric_limits<double>::infinity();

    for(map<StateId, Token*>::iterator iter = cur_toks.begin();
        iter != cur_toks.end(); ++iter)
    {
        if(best_cost > iter->second->cost)
        {
            best_cost = iter->second->cost;
            best_token = iter->second;
        }
    }

    vector<int> result_rev;
    Token* path = best_token;
    while(path != NULL)
    {
        if(path->arc.olabel != 0)
        {
            result_rev.push_back(path->arc.olabel);
        }
        path = path->prev;
    }

    vector<int> result;
    for(int i=result_rev.size()-1; i>=0; i--)
    {
        result.push_back(result_rev[i]);
    }

    return result;
}

void SimpleDecoder::AdvanceDecoding(P_Matrix feature, float acoustic_scale)
{
    while (num_frames_decoded < feature->rows)
    {
        // note: ProcessEmitting() increments num_frames_decoded_
        ClearToks(prev_toks);
        cur_toks.swap(prev_toks);
        ProcessEmitting(feature, acoustic_scale);
        ProcessNonemitting();
        PruneToks(m_beam, &cur_toks);
    }
}

void SimpleDecoder::ProcessEmitting(P_Matrix feature, float acoustic_scale)
{
    int32 frame = num_frames_decoded;
    // Processes emitting arcs for one frame.  Propagates from
    // prev_toks_ to cur_toks_.
    double cutoff = numeric_limits<float>::infinity();
    for(map<StateId, Token*>::iterator iter = prev_toks.begin();
       iter != prev_toks.end();
       ++iter)
    {
        StateId state = iter->first;
        Token *tok = iter->second;
        for(int i=0; i<m_fst->state[state].arcNum; i++)
        {
            P_Arc arc = &m_fst->state[state].arc[i];
            if(arc->ilabel != 0)
            {
                // propagate..
                float acoustic_cost = -acoustic_scale * LogLikelihood(feature, frame, arc->ilabel);
                double total_cost = tok->cost + arc->weight + acoustic_cost;

                if(total_cost >= cutoff)
                {
                    continue;
                }
                if(total_cost + m_beam < cutoff)
                {
                    cutoff = total_cost + m_beam;
                }

                Token *new_tok = newToken(arc, acoustic_cost, tok);
                map<StateId, Token*>::iterator find_iter = cur_toks.find(arc->nextstate);
                if(find_iter == cur_toks.end())
                {
                    cur_toks[arc->nextstate] = new_tok;
                }
                else
                {
                    if(find_iter->second->cost > new_tok->cost)
                    {
                        TokenDelete(find_iter->second);
                        find_iter->second = new_tok;
                    }
                    else
                    {
                        TokenDelete(new_tok);
                    }
               }
            }
        }
    }
    num_frames_decoded++;
}

void SimpleDecoder::ProcessNonemitting()
{
    // Processes nonemitting arcs for one frame.  Propagates within
    // cur_toks_.
    vector<StateId> queue;
    double infinity = std::numeric_limits<double>::infinity();
    double best_cost = infinity;

    for(map<StateId, Token*>::iterator iter = cur_toks.begin();
        iter != cur_toks.end();
        ++iter)
    {
        queue.push_back(iter->first);
        best_cost = min(best_cost, iter->second->cost);
    }
    double cutoff = best_cost + m_beam;

    while(!queue.empty())
    {
        StateId state = queue.back();
        queue.pop_back();
        Token *tok = cur_toks[state];
        for(int i=0; i<m_fst->state[state].arcNum; i++)
        {
            P_Arc arc = &m_fst->state[state].arc[i];

            if(arc->ilabel == 0)
            {   // propagate nonemitting only...
                const float acoustic_cost = 0.0;
                Token *new_tok = newToken(arc, acoustic_cost, tok);
                if(new_tok->cost > cutoff)
                {
                    TokenDelete(new_tok);
                }
                else
                {
                    map<StateId, Token*>::iterator find_iter = cur_toks.find(arc->nextstate);
                    if(find_iter == cur_toks.end())
                    {
                        cur_toks[arc->nextstate] = new_tok;
                        queue.push_back(arc->nextstate);
                    }
                    else
                    {
                        if(find_iter->second->cost > new_tok->cost)
                        {
                            TokenDelete(find_iter->second);
                            find_iter->second = new_tok;
                            queue.push_back(arc->nextstate);
                        }
                        else
                        {
                            TokenDelete(new_tok);
                        }
                    }
                }
            }
        }
    }
}

static const float kMinLogDiffFloat = logf(FLT_EPSILON);

float LogSumExp(vector<float> input, float prune)
{
    float max_elem = input[0];
    for(int i=1; i<input.size(); i++)
    {
        if(max_elem < input[i])
        {
            max_elem = input[i];
        }
    }
    float cutoff;
    cutoff = max_elem + kMinLogDiffFloat;
    if (prune > 0.0 && max_elem - prune > cutoff) // explicit pruning...
    {
        cutoff = max_elem - prune;
    }

    double sum_relto_max_elem = 0.0;

    for(int i = 0; i < input.size(); i++)
    {
        float f = input[i];
        if (f >= cutoff)
        {
            sum_relto_max_elem += expf(f - max_elem);
        }
    }
    return max_elem + logf(sum_relto_max_elem);
}

float SimpleDecoder::LogLikelihood(P_Matrix feature, int32 frame, int32 tid)
{
    int32 state = m_transmodel->TransitionIdToPdf(tid);

    vector<float> data;
    vector<float> data_squared;
    for(int i=0; i<feature->cols; i++)
    {
        float v = ReadMatrix(feature, frame, i);
        data.push_back(v);
        data_squared.push_back(v*v);
    }

    DiagGmm& pdf = m_amgmm->GetPdf(state);

    vector<float> loglikes;
    for(int i=0; i<pdf.gconsts.size(); i++)
    {
        loglikes.push_back(pdf.gconsts[i]);
    }

    for(int i=0; i<loglikes.size(); i++)
    {
        float sum = 0.0;
        for(int j=0; j<pdf.means_invvars.cols; j++)
        {
            sum += 1.0 * data[j] * ReadMatrix(&pdf.means_invvars, i, j);
            sum += -0.5 * data_squared[j] * ReadMatrix(&pdf.inv_vars, i, j);
        }
        loglikes[i] += sum;
    }

    float log_sum_exp_prune = -1;

    return LogSumExp(loglikes, log_sum_exp_prune);
}

void SimpleDecoder::ClearToks(map<StateId, Token*> &toks) {
    for(map<StateId, Token*>::iterator iter = toks.begin();
        iter != toks.end(); ++iter)
    {
        TokenDelete(iter->second);
    }
    toks.clear();
}

void SimpleDecoder::PruneToks(float beam, map<StateId, Token*> *toks)
{
    if(toks->empty())
    {
        printf("No tokens to prune.\n");
        return;
    }
    double best_cost = numeric_limits<double>::infinity();
    for(map<StateId, Token*>::iterator iter = toks->begin();
       iter != toks->end(); ++iter)
    {
        best_cost = min(best_cost, iter->second->cost);
    }

    vector<StateId> retained;
    double cutoff = best_cost + beam;
    for(map<StateId, Token*>::iterator iter = toks->begin();
       iter != toks->end(); ++iter)
    {
        if(iter->second->cost < cutoff)
        {
            retained.push_back(iter->first);
        }
        else
        {
            TokenDelete(iter->second);
        }
    }
    map<StateId, Token*> tmp;
    for (size_t i = 0; i < retained.size(); i++)
    {
        tmp[retained[i]] = (*toks)[retained[i]];
    }
    printf("Pruned to %lu toks.\n", retained.size());
    tmp.swap(*toks);
}