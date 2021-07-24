#include "common.h"

void ReadToken(FILE *fp, char* s)
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

void ReadIntegerVector(FILE *fp, vector<int32> *v)
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

void ReadBasicType(FILE *fp, int32 *t)
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

void ReadBasicType(FILE *fp, BaseFloat *t)
{
    uint8 size = 0;
    fread(&size, sizeof(size), 1, fp);

    if(size != sizeof(BaseFloat))
    {
        printf("float size error!\n");
        return;
    }

    fread(t, sizeof(*t), 1, fp);
}

void ReadFloatVectors(FILE *fp, vector<BaseFloat> *v)
{
    //TODO: Support other type, eg, double
    const char *my_token = "FV";
    char token[128];
    ReadToken(fp, token); //FV
    int32 size;
    ReadBasicType(fp, &size);
    v->resize(size);
    fread(v->data(), sizeof(BaseFloat), size, fp);
}

void ReadFloatMatrix(FILE *fp, P_Matrix m)
{
    const char *my_token = "FM";
    char token[128];
    ReadToken(fp, token); //FM

    int32 rows, cols;
    ReadBasicType(fp, &rows);
    ReadBasicType(fp, &cols);

    m->rows = rows;
    m->cols = cols;

    int32 skip = ((16 / sizeof(BaseFloat)) - cols % (16 / sizeof(BaseFloat))) % (16 / sizeof(BaseFloat));
    m->stride = cols + skip;

    int32 size = rows * cols;
    m->data.resize(size);
    fread(m->data.data(), sizeof(BaseFloat), size, fp);
}

BaseFloat ReadMatrix(P_Matrix m, int32 row, int32 col)
{
    return m->data[m->cols*row+col];
}