#include "fstreader.h"

const int32 fstMagicNumber = 2125659606;

bool FstHeader::Read(const char* fileName)
{
    printf("%s\n", fileName);
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

    ReadString(fsttype, fp);
    ReadString(arctype, fp);
    ReadInt(&version, sizeof(version), fp);
    ReadInt(&flags, sizeof(flags), fp);
    ReadInt(&properties, sizeof(properties), fp);
    ReadInt(&start, sizeof(start), fp);
    ReadInt(&numstates, sizeof(numstates), fp);
    ReadInt(&numarcs, sizeof(numarcs), fp);

    fclose(fp);
    return true;
}

void FstHeader::ReadString(char *buf, FILE *fp)
{
    uint32 len = 0;
    fread(&len, sizeof(len), 1, fp);
    buf = (char*)malloc(len+1);
    memset(buf, 0, len+1);
    fread(buf, len, 1, fp);
}

void FstHeader::ReadInt(void *buf, int bytes, FILE *fp)
{
    fread(buf, bytes, 1, fp);
}