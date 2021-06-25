#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

#define SAFE_FREE(x) if(x) {free(x); x=NULL;}

#endif