#include "fstreader.h"
#include "simple-decoder.h"

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        printf("arg error\n");
        return -1;
    }

    char* fstFileName = argv[1];

    float beam = 16.0;

    FstReader fstReader;
    fstReader.Read(fstFileName);

    SimpleDecoder decoder(&fstReader, beam);

    decoder.Decode();

    return 0;
}