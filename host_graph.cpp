
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <pthread.h>

#include "host_graph_sw.h"


#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/random.hpp>

using namespace boost::multiprecision;
using namespace boost::random;



#define RNG_RANGE_MAX           (30)


//
// Generate the numbers:
//

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
} thread_item;

#define MAX_NUN_THREADS 16
thread_item   threads[MAX_NUN_THREADS];

#define  MAX_NUM_ESTIMATOR  (1000 * 1000 * 10)

#define STATUS_FIRST_EDGE       (0)
#define STATUS_SECOND_EDGE      (1)
#define STATUS_CLOSE            (2)
#define STATUS_END              (3)


estimator local_estimator[MAX_NUM_ESTIMATOR];

uniform_int_distribution<mpz_int> ui(0, mpz_int(1) << RNG_RANGE_MAX);

const double ratio = 1.0 / (double(1 << RNG_RANGE_MAX));

typedef mt11213b  prng;

int reservoir_sampling(int n, prng &lmt)
{
    bernoulli_distribution<double> p(1.0 / (n));
    bool res = p(lmt);
    if (res)
    {
        return 1;
    }
    return 0;
}



int approximation(estimator *p_est, Graph* gptr, CSR* csr)
{
    prng mt;
    mt.seed(static_cast<unsigned int>(std::time(0)) + (unsigned long)(p_est));
    memset(p_est, 0, sizeof(estimator));
    int edgeNum   = csr ->edgeNum;

    sample_edge * p_first    = &p_est->first_edge;
    sample_edge * p_second    = &p_est->second_edge;
    p_est->status = 0;

    for (int i = 0; i < edgeNum; i++)
    {
        if (reservoir_sampling(i + 1, mt) )
        {
            p_est->first_edge.node[0] = gptr->data[i][0];
            p_est->first_edge.node[1] = gptr->data[i][1];
            p_first->update_counter ++;
            p_est->neighbor_counter = 0;
            p_est->status = 1;
        }
        else
        {
            int ln = gptr->data[i][0];
            int rn = gptr->data[i][1];


            int adjacent_flag = 0;
#if 0
            if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
            {
                p_est->neighbor_flag = 1;
                adjacent_flag = 1;
            }
#endif
#if 1
            if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
            {
                p_est->neighbor_flag = 0;
                adjacent_flag = 1;
            }
#endif
            if (adjacent_flag == 1)
            {
                p_est->neighbor_counter ++;
                if (reservoir_sampling(p_est->neighbor_counter, mt))
                {
                    p_second->node[0] =  (p_first->node[1] == ln) ? ln : rn;
                    p_second->node[1] =  (p_first->node[1] == ln) ? rn : ln;
                    p_second->update_counter ++;
                    p_est->status = 2;
                }

                if (p_est->status == 2)
                {
#if 1
                    int node = p_first->node[p_est->neighbor_flag ] ;
                    //printf("node %d -> %d\n", node, p_second_edge->node[1]);

                    for (int j = csr->rpao[node]; j < csr->rpao[node + 1]; j++)
                    {
                        //printf("%d\n", csr->ciao[j]);
                        if (csr->ciao[j] == p_second->node[1])
                        {
                            p_est->expecation = p_est->neighbor_counter * i;
                            p_est->status = 3;
                            break;
                        }
                    }
#else
                    if (((ln == p_first->node[p_est->neighbor_flag]) && (rn == p_second->node[1]))
                            || ((rn == p_first->node[p_est->neighbor_flag]) && (ln == p_second->node[1])))

                    {

                        p_est->expecation = p_est->neighbor_counter  * i;
                        p_est->status = 3;
                        return 1;
                    }
#endif
                }

            }
        }
    }

    return  0;
}

void *approximation_thread(void *argv)
{

    //DEBUG_PRINTF("here\n");
    thread_item *p_data = (thread_item *)argv;
    approximation(p_data->est, p_data->gptr, p_data->csr);
    return 0;
}

using namespace std;

graphInfo graphDataInfo;

#define EST_NUM                 (10)

int main(int argc, char **argv) {

    std::string gName;
    if (argc > 2)
    {
        gName = argv[2];
    }
    else
    {
        gName = "rmat-19-32";
    }

    int est_num = (argc < 4) ? EST_NUM : std::stoi(argv[3]);
    std::string mode = "normal";
#if 0
    Graph* gptr = createGraph(gName, mode);
    CSR* csr    = new CSR(*gptr);
    csr->save2File(gName);
    free(gptr);
    return 0;
#endif

    DEBUG_PRINTF("start main\n");

    std::cout << std::hex << std::showbase;
    std::cout << (mpz_int(1) << RNG_RANGE_MAX) << std::endl;

    Graph* gptr = createGraph(gName, mode);
    CSR* csr    = new CSR(*gptr);

    int edgeNum   = csr ->edgeNum;


//ARRAY_SIZE(local_estimator)
    int est_num_map [] = { 16, 100, 1000, 2000,
                           5000, 10000, 20000, 50000,
                           100000, 200000, 500000, 1000 * 1000,
                           1000 * 1000 * 2, 1000 * 1000 * 5
                         };
    //ARRAY_SIZE(est_num_map)
    for (int k = 0; k < ARRAY_SIZE(est_num_map); k ++)
    {
        est_num = est_num_map[k];
        for (int repeat = 0; repeat < 1; repeat ++)
        {

            for (int j = 0; j < ARRAY_SIZE(threads); j++)
            {
                threads[j].expecation = 0.0;
                threads[j].counter = 0;
                threads[j].run = 0;
            }

            for (int i = 0 ; i < est_num / 16; i ++)
            {
                for (int j = 0; j < ARRAY_SIZE(threads); j++)
                {
                    threads[j].est  = &local_estimator[j];
                    threads[j].gptr = gptr;
                    threads[j].csr  = csr;
                    pthread_create(&threads[j].pid, NULL, approximation_thread, &threads[j]);
                }

                for (int j = 0; j < ARRAY_SIZE(threads); j++)
                {
                    pthread_join(threads[j].pid, NULL);
                    threads[j].expecation += local_estimator[j].expecation;
                    if (local_estimator[j].status == 3)
                    {
                        threads[j].counter++;
                    }
                    threads[j].run ++;
                    estimator * p_est    = &local_estimator[j];
                    sample_edge * p_first     = &p_est->first_edge;
                    sample_edge * p_second = &p_est->second_edge;
#if 0
                    DEBUG_PRINTF("%d nc %d (%d %d) [%d %d] [%d %d]\n",
                                 p_est->status ,
                                 p_est->neighbor_counter,
                                 p_first->update_counter,
                                 p_second->update_counter,
                                 p_first->node[0],
                                 p_first->node[1],
                                 p_second->node[0],
                                 p_second->node[1]);
#endif
                }

                //estimator * p_est    = &local_estimator[i];

                //success_counter += approximation(p_est, gptr, csr);
            }
            double result = 0;
            int total_est_num = 0;
            int success_counter = 0;

            for (int j = 0; j < ARRAY_SIZE(threads); j++)
            {
                result += threads[j].expecation;
                total_est_num += threads[j].run;
                success_counter += threads[j].counter;

            }
#if 0
            for (int i = 0; i < est_num; i++)
            {

                estimator * p_est    = &local_estimator[i];
                sample_edge * p_edge    = &p_est->first_edge;
                sample_edge * p_second_edge = &p_est->second_edge;

                DEBUG_PRINTF("id %d nc %d %d status %d [%d %d] [%d %d]\n",
                             p_edge->id,
                             p_est->neighbor_counter,
                             p_est->neighbor_id,
                             p_second_edge->id,
                             p_edge->node[0],
                             p_edge->node[1],
                             p_second_edge->node[0],
                             p_second_edge->node[1]);

                result     += p_est->expecation;
            }
#endif
            DEBUG_PRINTF("result %lf @ %d with %d,total %d \n", (double(result ) / total_est_num),
                         repeat, total_est_num,
                         success_counter);
        }
    }
//    for (unsigned i = 0; i < 10; ++i)
//        std::cout << ui(mt) << std::endl;

    free(gptr);
    free(csr);

    return 0;
}

