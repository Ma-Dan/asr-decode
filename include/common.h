#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>

using namespace std;

typedef unsigned char uint8;
typedef signed char int8;
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

#define SAFE_FREE(x) if(x) {free(x); x=NULL;}

#ifndef M_LOG_2PI
#define M_LOG_2PI 1.8378770664093454835606594728112
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON 1.19209290e-7f
#endif

typedef struct tagMatrix
{
    int32 cols;
    int32 rows;
    int32 stride;
    vector<float> data;
} Matrix, *P_Matrix;

void ReadToken(FILE *fp, char* s);
void ReadIntegerVector(FILE *fp, vector<int32> *v);
void ReadBasicType(FILE *fp, int32 *t);
void ReadBasicType(FILE *fp, float *t);
void ReadFloatVectors(FILE *fp, vector<float> *v);
void ReadFloatMatrix(FILE *fp, P_Matrix m);
float ReadMatrix(P_Matrix m, int32 row, int32 col);

#endif