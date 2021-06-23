#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>

#define SAFE_FREE(x) if(x) {free(x); x=NULL;}

#endif