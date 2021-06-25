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
    return arc->weight1;// + arc->weight2;
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

SimpleDecoder::SimpleDecoder(FstReader *fst, float beam)
{
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

bool SimpleDecoder::Decode()
{
    InitDecoding();
    AdvanceDecoding();
    return (!cur_toks.empty());
}

void SimpleDecoder::AdvanceDecoding()
{
    while (num_frames_decoded < 668)
    {
        // note: ProcessEmitting() increments num_frames_decoded_
        ClearToks(prev_toks);
        cur_toks.swap(prev_toks);
        ProcessEmitting();
        ProcessNonemitting();
        PruneToks(m_beam, &cur_toks);
  }
}

void SimpleDecoder::ProcessEmitting()
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
                // TODO: Add AM cost here
                float acoustic_cost = 0;
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