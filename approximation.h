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


#define MAX_NUM_ESTIMATOR   (1000 * 1000 * 10)
#define MAX_NUN_THREADS     (10)


typedef struct
{
    unsigned int id;
    int node[2];
    int temp_deg[2];
    unsigned int p;
    unsigned int update_counter;
} sample_edge;


typedef uint64_t   uint64_rng_t;
typedef uint32_t   uint32_rng_t;


struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_rng_t state;             // RNG state.  All values are possible.
    uint64_rng_t inc;               // Controls which RNG sequence (stream) is
                                   // selected. Must *always* be odd.
};

typedef struct pcg_state_setseq_64 pcg32_random_t;

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
    unsigned int sub_est_num;
    unsigned int temp_neighbor_counter;
    unsigned int success;

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


void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq);


typedef mt11213b  prng;

int approximation_motifs_scheme_1(estimator *p_est, Graph* gptr, int edgeNum, int id);
int approximation_motifs_scheme_2(estimator *p_est, Graph* gptr, int edgeNum, int id);
int approximation_motifs_scheme_3(estimator *g_est, Graph* gptr, int edgeNum , int est_id);
int approximation_motifs_scheme_4(estimator *g_est, Graph* gptr, int edgeNum , int est_id);

int triangle_count(Graph* gptr, CSR* csr, int num);

int approximation_triangle_scheme_1(estimator *p_est, Graph* gptr, int edgeNum, int id);
int approximation_triangle_scheme_2(estimator *p_est, Graph* gptr, int edgeNum, int id);
int approximation_triangle_scheme_3(estimator *g_est, Graph* gptr, int edgeNum , int est_id);
int approximation_triangle_scheme_4(estimator *g_est, Graph* gptr, int edgeNum , int est_id);

extern prng  mt;

extern pthread_mutex_t lock;


int pcg_reservoir_sampling(int n, pcg32_random_t* rng);
int reservoir_sampling(int n, prng &lmt);


int rng_pcg_res_test(int edgeNum, int num);
int rng_res_test(int edgeNum, int num);
int rng_test(int edgeNum, int num);




#endif
