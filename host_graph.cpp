
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <pthread.h>

#include "host_graph_sw.h"

#include "approximation.h"


#define APPROXIMATE_FUNCTION         approximation_triangle_scheme_1

//
// Generate the numbers:
//


#define MAX_NUN_THREADS 16
thread_item   threads[MAX_NUN_THREADS];

#define  MAX_NUM_ESTIMATOR  (1000 * 1000 * 10)

#define STATUS_FIRST_EDGE       (0)
#define STATUS_SECOND_EDGE      (1)
#define STATUS_CLOSE            (2)
#define STATUS_END              (3)


estimator local_estimator[MAX_NUN_THREADS];


void *approximation_thread(void *argv)
{

    //DEBUG_PRINTF("here\n");
    thread_item *p_data = (thread_item *)argv;
    APPROXIMATE_FUNCTION(p_data->est, p_data->gptr, p_data->csr, p_data->est_id);
    return 0;
}


int rng_test(CSR* csr, int num)
{
    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0))  + counter);

    counter ++;
    int edgeNum   = csr ->edgeNum;
    uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);
    for (int i = 0; i < num; i ++)
    {
        int loc = static_cast<int>(ui(mt));
        DEBUG_PRINTF("rng %d %d \n", loc , edgeNum);
    }
    return 0;
}

int rng_res_test(CSR* csr, int num)
{
    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0))  + counter);

    counter ++;
    int edgeNum   = csr ->edgeNum;
    uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);
    for (int i = 0; i < num; i ++)
    {
        int loc = 0;
        for (int k = 0; k < edgeNum; k ++)
        {
            if (reservoir_sampling(k + 1, mt) )
            {
                loc = k;
            }
        }

        DEBUG_PRINTF("rng %d %d \n", loc , edgeNum);
    }
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
    int printf_flag = std::stoi(argv[1]);
    std::string mode = "normal";

    DEBUG_PRINTF("start main\n");

    std::cout << std::hex << std::showbase;

    Graph* gptr = createGraph(gName, mode);
    CSR* csr    = new CSR(*gptr);

    int edgeNum   = csr ->edgeNum;

    //rng_res_test(csr, 100000);
    //return 0;

   //1465313402


    int est_num_map [] = { 16, 100, 1000, 2000,
                           5000, 10000, 20000, 50000,
                           100000, 200000, 500000, 1000 * 1000,
                           1000 * 1000 * 2, 1000 * 1000 * 5
                         };
#if 1
    const int total_k = ARRAY_SIZE(est_num_map) - 3;
    for (int k = 0; k < total_k ; k ++)
    {
        est_num = est_num_map[k];
        triangle_count(gptr, csr, est_num);
    }

    free(gptr);
    free(csr);
    return 0;
#else
    //ARRAY_SIZE(est_num_map)
    int est_id = 0;
    const int total_k = ARRAY_SIZE(est_num_map);
    for (int k = 0; k < total_k ; k ++)
    {
        est_num = est_num_map[k];
        for (int repeat = 0; repeat < 1; repeat ++)
        {

            for (int j = 0; j < ARRAY_SIZE(threads); j++)
            {
                threads[j].expecation = 0.0;
                threads[j].counter = 0;
                threads[j].run = 0;
                est_id ++;
                threads[j].est_id = est_id;
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

                    if (local_estimator[j].status == 3)
                    {
                        threads[j].expecation += local_estimator[j].expecation;
                        threads[j].counter++;
                    }
                    threads[j].run ++;
                    estimator * p_est    = &local_estimator[j];
                    sample_edge * p_first     = &p_est->first_edge;
                    sample_edge * p_second = &p_est->second_edge;
                    if  (printf_flag == 1)
                    {
                        DEBUG_PRINTF("%d nc %d (%d %d) %d-[%d %d] %d-[%d %d] --> %f\n",
                                     p_est->status ,
                                     p_est->neighbor_counter,
                                     p_first->update_counter,
                                     p_second->update_counter,
                                     p_first->id,
                                     p_first->node[0],
                                     p_first->node[1],
                                     p_second->id,
                                     p_second->node[0],
                                     p_second->node[1],
                                     p_est->expecation);
                    }
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
            if (1)
            {
                DEBUG_PRINTF("result %lf @ %d with %d,total %d \n", (double(result ) / total_est_num),
                             repeat, total_est_num,
                             success_counter);

            }
        }

    }
//    for (unsigned i = 0; i < 10; ++i)
//        std::cout << ui(mt) << std::endl;

    free(gptr);
    free(csr);

    return 0;
#endif
}

