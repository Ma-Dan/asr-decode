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

#endif