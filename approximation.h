#ifndef __APPROXIMATION_H__
#define __APPROXIMATION_H__

//#include "host_graph_sw.h"
#include "graph.h"
#include "common.h"
#include "pthread.h"

#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/random.hpp>

using namespace boost::multiprecision;
using namespace boost::random;


#define  MAX_NUM_ESTIMATOR  (1000 * 1000 * 10)


typedef struct
{
    unsigned int id;
    int node[2];
    unsigned int p;
    unsigned int update_counter;
} sample_edge;



typedef struct {
    double expecation;
    sample_edge first_edge;
    sample_edge second_edge;
    unsigned int tmp_neighbor_counter;
    unsigned int tmp_rng;
    unsigned int neighbor_counter;
    unsigned int neighbor_id;
    unsigned char neighbor_flag;
    unsigned int status;

} estimator;

typedef struct
{
    pthread_t pid;
    estimator *est;
    Graph* gptr;
    CSR* csr;
    double expecation;
    unsigned int run;
    unsigned int counter;
    unsigned int est_id;
    int edgeNum;
} thread_item;


typedef mt11213b  prng;

int approximation_motifs_scheme_1(estimator *p_est, Graph* gptr, int edgeNum, int id);
int approximation_motifs_scheme_2(estimator *p_est, Graph* gptr, int edgeNum, int id);


int triangle_count(Graph* gptr, CSR* csr, int num);

int approximation_triangle_scheme_1(estimator *p_est, Graph* gptr, int edgeNum, int id);
int approximation_triangle_scheme_2(estimator *p_est, Graph* gptr, int edgeNum, int id);

extern prng  mt;

extern pthread_mutex_t lock;

inline int reservoir_sampling(int n, prng &lmt)
{
    boost::random::bernoulli_distribution<double> p(1.0 / (n));
    pthread_mutex_lock(&lock);
    bool res = p(mt);
    pthread_mutex_unlock(&lock);
    //DEBUG_PRINTF("n %d, res %d\n", n, res);
    if (res)
    {
        return 1;
    }
    return 0;
}





#endif
