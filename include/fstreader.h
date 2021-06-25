#ifndef FST_READER
#define FST_READER

#include "common.h"

typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

class FstHeader {
    public:
        bool Read(const char* fileName);
        ~FstHeader();

        char* fsttype;      // E.g. "vector".
        char* arctype;      // E.g. "standard".
        int32 version;      // Type version number.
        int32 flags;        // File format bits.
        uint64 properties;  // FST property bits.
        int64 start;        // Start state.
        int64 numstates;    // # of states.
        int64 numarcs;      // # of arcs.
    private:
        void ReadString(char **buf, FILE *fp);
        void ReadInt(void *buf, int bytes, FILE *fp);
};

typedef struct tagArc
{
    int ilabel;
    int olabel;
    float weight;
    int nextstate;
} Arc, *P_Arc;

typedef struct tagState
{
    float weight;
    int field1;
    int arcNum;
    int field3;
    int field4;
    P_Arc arc;
} State, *P_State;

class FstReader {
    public:
        bool Read(const char* fileName);
        int Start();
        ~FstReader();
    //private:
        FstHeader hdr;
        P_State state;
        P_Arc arc;
};

#endif