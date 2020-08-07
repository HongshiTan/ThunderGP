
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <pthread.h>

#include "approximation.h"

#define STR_APPROXIMATE_FUNCTION        STRINGIFY_MACRO(APPROXIMATE_FUNCTION)

#define  MAX_NUM_ESTIMATOR  (1000 * 1000 * 10)

estimator local_estimator[MAX_NUN_THREADS];
thread_item   threads[MAX_NUN_THREADS];


void *approximation_thread(void *argv)
{
    thread_item *p_data = (thread_item *)argv;
    APPROXIMATE_FUNCTION(p_data->est, p_data->gptr, p_data->edgeNum, p_data->est_id);
    return 0;
}


using namespace std;


prng mt;
pthread_mutex_t lock;

int main(int argc, char **argv) {

    mt.seed(static_cast<unsigned int>(std::time(0)));

    std::string gName;
    if (argc > 2)
    {
        gName = argv[2];
    }
    else
    {
        gName = "rmat-19-32";
    }

    int est_num = (argc < 4) ? 1000 : std::stoi(argv[3]);
    int printf_flag = std::stoi(argv[1]);
    std::string mode = "normal";

    DEBUG_PRINTF("start main: %s \n", STR_APPROXIMATE_FUNCTION);

    std::cout << std::hex << std::showbase;

    Graph* gptr = createGraph(gName, mode);
    DEBUG_PRINTF("edge num :%d \n", gptr->edgeNum);

    int edgeNum   = gptr ->edgeNum;


    if (pthread_mutex_init(&lock, NULL) != 0)
    {
        printf("\n mutex init failed\n");
        return 1;
    }

    int est_num_map [] = {  10,       20,       50, 
                            100,      200,      500, 
                            1000,     2000,     5000, 
                            10000,    20000,    50000,
                            100000,   200000,   500000,
                            1000000,  2000000,  5000000,
                            10000000, 20000000, 50000000,
                         };
    int est_id = 0;

#ifndef MAX_RUN_STEPS 
    const int total_k = ARRAY_SIZE(est_num_map);
#else
    const int total_k = MAX_RUN_STEPS;
#endif

    for (int k = 0; k < total_k ; k ++)
    {
        est_num = est_num_map[k];
        boost::random::uniform_int_distribution<mpz_int> iui(0, (edgeNum - 1) );


        for (int j = 0; j < ARRAY_SIZE(threads); j++)
        {
            threads[j].expecation = 0.0;
            threads[j].counter = 0;
            threads[j].run = 0;
            est_id ++;
            threads[j].est_id = est_id + static_cast<int>(iui(mt));
        }

        for (int i = 0 ; i < est_num / ARRAY_SIZE(threads); i ++)
        {
            //DEBUG_PRINTF("%d i \n",i);
            for (int j = 0; j < ARRAY_SIZE(threads); j++)
            {
                local_estimator[j].sub_est_num = SUB_EST;
                threads[j].est  = &local_estimator[j];
                threads[j].gptr = gptr;
                threads[j].edgeNum  = gptr->edgeNum;
                pthread_create(&threads[j].pid, NULL, approximation_thread, &threads[j]);
            }

            for (int j = 0; j < ARRAY_SIZE(threads); j++)
            {
                pthread_join(threads[j].pid, NULL);

                if (local_estimator[j].status == 3)
                {
                    threads[j].expecation += local_estimator[j].expecation;
                    threads[j].counter    += local_estimator[j].success;
                }
                threads[j].run ++;
                estimator * p_est      = &local_estimator[j];
                sample_edge * p_first  = &p_est->first_edge;
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
            DEBUG_PRINTF("result %lf - %lf  with %d, success %d ratio %lf %lf\n", 
                        (double(result) * ((double)edgeNum / (total_est_num * SUB_EST))),
                         double(result),
                         total_est_num * SUB_EST,
                         success_counter, ((double)success_counter / (total_est_num * SUB_EST)), ((double)result / (double)success_counter));

        }


    }

    free(gptr);
    return 0;
}

