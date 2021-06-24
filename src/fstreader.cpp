#include "fstreader.h"

const int32 fstMagicNumber = 2125659606;

bool FstHeader::Read(const char* fileName)
{
    FILE *fp;
    fp = fopen(fileName, "rb");

    if(fp == NULL)
    {
        printf("Error opening fst file\n");
        return false;
    }

    int32 magic_number = 0;
    ReadInt(&magic_number, sizeof(magic_number), fp);
    if (magic_number != fstMagicNumber) {
        printf("FstHeader::Read: Bad FST header\n");
        return false;
    }

    ReadString(&fsttype, fp);
    ReadString(&arctype, fp);
    ReadInt(&version, sizeof(version), fp);
    ReadInt(&flags, sizeof(flags), fp);
    ReadInt(&properties, sizeof(properties), fp);
    ReadInt(&start, sizeof(start), fp);
    ReadInt(&numstates, sizeof(numstates), fp);
    ReadInt(&numarcs, sizeof(numarcs), fp);

    fclose(fp);
    return true;
}

void FstHeader::ReadString(char **buf, FILE *fp)
{
    uint32 len = 0;
    fread(&len, sizeof(len), 1, fp);
    *buf = (char*)malloc(len+1);
    memset(*buf, 0, len+1);
    fread(*buf, len, 1, fp);
}

void FstHeader::ReadInt(void *buf, int bytes, FILE *fp)
{
    fread(buf, bytes, 1, fp);
}

FstHeader::~FstHeader()
{
    SAFE_FREE(fsttype);
    SAFE_FREE(arctype);
}


bool FstReader::Read(const char* fileName)
{
    if(!hdr.Read(fileName))
    {
        return false;
    }

    FILE *fp;
    fp = fopen(fileName, "rb");

    if(fp == NULL)
    {
        printf("Error opening fst file\n");
        return false;
    }

    //65 bytes header
    fseek(fp, 65, SEEK_SET);
    //Check the type of Arc

    //Read the FST
    //20 bytes per state
    state = (P_State)malloc(hdr.numstates * sizeof(State));
    for(int64 i=0; i<hdr.numstates; i++)
    {
        fread(&state[i], 20, 1, fp);
    }

    //16 bytes per arc
    arc = (P_Arc)malloc(hdr.numarcs * sizeof(Arc));
    fread(arc, sizeof(Arc), hdr.numarcs, fp);

    //Assign arcs to states
    int64 offset = 0;
    for(int64 i=0; i<hdr.numstates; i++)
    {
        state[i].arc = &arc[offset];
        offset += state[i].arcNum;
    }

    fclose(fp);
    return true;
}

FstReader::~FstReader()
{
    SAFE_FREE(state);
    SAFE_FREE(arc);
}