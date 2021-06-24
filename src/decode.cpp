#include "fstreader.h"

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        printf("arg error\n");
        return -1;
    }

    char* fstFileName = argv[1];

    FstReader fstReader;
    fstReader.Read(fstFileName);
    return 0;
}