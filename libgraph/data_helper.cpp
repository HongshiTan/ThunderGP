#include <stdio.h>

#include "graph.h"
#include "common.h"
//#include "host_graph_sw.h"

static int startIdx = 0;



int localGetStartIndex(const std::string &name)
{
    int startVertexIdx = 1;
    if (name == "youtube")    startVertexIdx = 320872;
    if (name == "lj1")        startVertexIdx = 3928512;
    if (name == "pokec")      startVertexIdx = 182045;
    if (name == "rmat-19-32") startVertexIdx = 104802;
    if (name == "rmat-21-32") startVertexIdx = 365723;

    return startVertexIdx;
}


Graph* createGraph(const std::string &gName, const std::string &mode) {
    Graph* gptr;
    std::string dir;
    if (mode == "normal") dir = "/graph_data/";

    else {
        std::cout << "unknown execution environment." << std::endl;
        exit(0);
    }

    gptr = new Graph(gName, 0);
    startIdx = localGetStartIndex(gName);
    return gptr;
}

int getStartIndex(void)
{
    return startIdx;
}



double getCurrentTimestamp(void) {
    timespec a;
    clock_gettime(CLOCK_MONOTONIC, &a);
    return (double(a.tv_nsec) * 1.0e-9) + double(a.tv_sec);
}

