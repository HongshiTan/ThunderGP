#ifndef __COMMON_H__
#define __COMMON_H__


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#include "config.h"

#define STRINGIFY_MACRO(x)              STR(x)
#define STR(x)                          #x
#define DUMP_CAPACITY(a)                DEBUG_PRINTF("mem size dump %s %u\n",#a, (uint32_t)a.capacity());

#if 1

#define DEBUG_PRINTF(fmt,...)   printf(fmt,##__VA_ARGS__); fflush(stdout);

#else

#define DEBUG_PRINTF(fmt,...)   ;

#endif


#ifndef ARRAY_SIZE
#define ARRAY_SIZE(arr) ((int)(sizeof(arr)/sizeof((arr)[0])))
#endif


#define SIZE_ALIGNMENT(in,align)     ((((((unsigned int)(in + 1))/(align)) + 1) * (align )) - 1) 

typedef struct
{
    int vertexNum;
    int compressedVertexNum;
    int edgeNum;
    int blkNum;
} graphInfo;

Graph* createGraph(const std::string &gName, const std::string &mode);

#endif /* __COMMON_H__ */
