#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

class FstHeader {
    public:
        bool Read(const char* fileName);
    private:
        char* fsttype;      // E.g. "vector".
        char* arctype;      // E.g. "standard".
        int32 version;      // Type version number.
        int32 flags;        // File format bits.
        uint64 properties;  // FST property bits.
        int64 start;        // Start state.
        int64 numstates;    // # of states.
        int64 numarcs;      // # of arcs.

        void ReadString(char *buf, FILE *fp);
        void ReadInt(void *buf, int bytes, FILE *fp);
};