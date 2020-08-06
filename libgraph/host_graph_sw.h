#ifndef __HOST_GRAPH_SW_H__
#define __HOST_GRAPH_SW_H__

#include <stdint.h>

#include "common.h"

#include "graph.h"

#include "he_mem.h"
#include "he_mem_id.h"



#define MAX_PARTITIONS_NUM      (128)

#define DEFAULT_KERNEL_ID       (0)


#define checkStatus(str) {                          \
    if (status != 0 || status != CL_SUCCESS) {      \
        printf("Error code: %d\n", status);         \
        printf("Error: %s\n", str);                 \
        acceleratorDeinit();                        \
        exit(-1);                                   \
    }                                               \
}



typedef struct
{
    double fpgaExeTime;
    double effic;
    double compress;
    double degree;
} profileLog;

typedef struct
{
    unsigned int  listStart;
    unsigned int  listEnd;
    unsigned int  dstVertexStart;
    unsigned int  dstVertexEnd;
    unsigned int  srcVertexStart;
    unsigned int  srcVertexEnd;
    float         scatterCacheRatio;
    float         compressRatio;
    he_mem_t      edge;
    he_mem_t      edgeMap;
    he_mem_t      edgeProp;
    he_mem_t      tmpProp;
    unsigned int  mapedTotalIndex;
    unsigned int  cuIndex;
    profileLog log;
} subPartitionDescriptor;

typedef struct
{
    subPartitionDescriptor  *sub[SUB_PARTITION_NUM];
    int                     finalOrder[SUB_PARTITION_NUM];
    unsigned int            totalEdge;
    unsigned int            subPartitionSize;
    cl_event                syncEvent[SUB_PARTITION_NUM];
    cl_event                applyEvent;
    double                  applyExeTime;
} partitionDescriptor;

typedef struct
{
    const char* name;
    int partition_mem_attr;
    int prop_id;
    int output_id;
    cl_kernel kernel;
    he_mem_t  prop[2];
    he_mem_t  tmpProp;
} gatherScatterDescriptor;

typedef struct
{
    const char* name;
    cl_kernel kernel;
} applyDescriptor;

typedef struct
{
    CSR* csr;

    subPartitionDescriptor subPartitions[MAX_PARTITIONS_NUM * SUB_PARTITION_NUM];

    partitionDescriptor partitions[MAX_PARTITIONS_NUM];

    gatherScatterDescriptor * gsKernel[SUB_PARTITION_NUM];

    applyDescriptor * applyKernel;

    cl_command_queue gsOps[SUB_PARTITION_NUM];

    cl_command_queue applyOps;

    cl_program program;

    cl_platform_id platform;

    cl_device_id device;

    cl_context context;
    
} graphAccelerator;


subPartitionDescriptor * getSubPartition(int partID);

partitionDescriptor * getPartition(int partID);



gatherScatterDescriptor * getGatherScatter(int kernelID);

applyDescriptor * getApply(void);

void kernelInit(graphAccelerator * acc);

void setGsKernel(int partId, int superStep, graphInfo *info);

void setApplyKernel(int partId, int superStep, graphInfo *info);



graphAccelerator * getAccelerator(void);

int acceleratorInit(const char * name, char *file_name);

int acceleratorDataPrepare(const std::string &gName, const std::string &mode, graphInfo *info);

int acceleratorDataPreprocess(graphInfo *info);

int acceleratorSuperStep(int superStep, graphInfo *info);

int accelratorProfile (int superStep, int runCounter, graphInfo *info, double exeTime);

int acceleratorDeinit(void);



void partitionGatherScatterCModel(
    cl_context              &context,
    cl_device_id            &device,
    int                     superStep,
    int                     cuIndex,
    subPartitionDescriptor  *subPartitions
);

void partitionApplyCModel(
    cl_context              &context,
    cl_device_id            &device,
    int                     superStep,
    int                     partId,
    unsigned int            applyArg
);

unsigned int dataPrepareGetArg(graphInfo *info);

int dataPrepareProperty(graphInfo *info);

int getStartIndex(void);




double getCurrentTimestamp(void);






#endif /* __HOST_GRAPH_SW_H__ */
