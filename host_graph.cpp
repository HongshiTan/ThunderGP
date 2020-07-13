
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>

#include "host_graph_sw.h"


#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/random.hpp>

using namespace boost::multiprecision;
using namespace boost::random;



#define RNG_RANGE_MAX           (30)

mt19937 mt;
uniform_int_distribution<mpz_int> ui(0, mpz_int(1) << RNG_RANGE_MAX);
//
// Generate the numbers:
//

typedef struct
{
    unsigned int id;
    int node[2];
    unsigned int p;
    unsigned int neighbor_counter;
} sample_edge;



typedef struct {
    double expecation;
    sample_edge first_edge;
    sample_edge second_edge;
    unsigned int neighbor_id;
    unsigned char neighbor_flag;
    unsigned int status;

} estimator;

#define  MAX_NUM_ESTIMATOR  (1000 * 1000 * 10)

#define STATUS_FIRST_EDGE       (0)
#define STATUS_SECOND_EDGE      (1)
#define STATUS_CLOSE            (2)
#define STATUS_END              (3)


estimator local_estimator[MAX_NUM_ESTIMATOR];


const double ratio = 1.0 / (double(1 << RNG_RANGE_MAX));


int reservoir_sampling(int n)
{
    int rn = static_cast<int>(ui(mt));
    float p = float(rn) * ratio;
    float p0 = 1.0f / (n + 1);
    if (p < p0)
    {
        return 1;
    }
    return 0;
}

using namespace std;

graphInfo graphDataInfo;

#define EST_NUM                 (10)

int main(int argc, char **argv) {

    char * xcl_file = NULL;
    if (argc > 1)
    {
        xcl_file = argv[1];
    }

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
    mt.seed(static_cast<unsigned int>(std::time(0)));

    std::cout << std::hex << std::showbase;
    std::cout << (mpz_int(1) << RNG_RANGE_MAX) << std::endl;

    Graph* gptr = createGraph(gName, mode);
    CSR* csr    = new CSR(*gptr);

    int vertexNum = csr ->vertexNum;
    int edgeNum   = csr ->edgeNum;


//ARRAY_SIZE(local_estimator)
    int est_num_map [] = { 10, 100, 1000, 2000,
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
            for (int i = 0 ; i < est_num; i ++)
            {
                estimator  * p_est =  &local_estimator[i];
                sample_edge * p_edge = &p_est->first_edge;
                int rn = static_cast<int>(ui(mt));
                int id = int((double (rn)) * ratio  * edgeNum);
                //printf("%d@%d %x id\n", id, edgeNum, rn);
                p_edge->id = id;
                p_edge->node[0] =  gptr->data[id][0];
                p_edge->node[1] =  gptr->data[id][1];
                p_edge->p = edgeNum;
                p_est->status = STATUS_SECOND_EDGE;
            }
            int success_counter = 0;

            for (int i = 0; i < est_num; i++)
            {
                estimator * p_est       = &local_estimator[i];
                sample_edge * p_edge    = &p_est->first_edge;
                sample_edge * p_second_edge = &p_est->second_edge;

                p_edge->neighbor_counter = 0;
                for (int j = p_edge->id + 1; j < edgeNum; j++)
                {
                    int ln =  gptr->data[j][0];
                    int rn =  gptr->data[j][1];
                    if ((ln == p_edge->node[1]) || (rn == p_edge->node[1]))
                    {
                        p_edge->neighbor_counter ++;
                        p_est->neighbor_flag = 0;
                    }
#if 0
                    else if ((ln == p_edge->node[0]) || (rn == p_edge->node[0]))
                    {
                        p_edge->neighbor_counter ++;
                        p_est->neighbor_flag = 1;
                    }
#endif
                }
                int rng = static_cast<int>(ui(mt));
                p_est->neighbor_id = int((double (rng)) * ratio  * p_edge->neighbor_counter);
                while (1)
                {
                    if(p_est->status== STATUS_END)
                    {
                        break;
                    }
                    switch (p_est->status)
                    {
                    case STATUS_SECOND_EDGE:
                    {

                        unsigned int local_neighbor_counter = 0;
                        p_est->status = STATUS_END;
                        for (int j = p_edge->id + 1; j < edgeNum; j++)
                        {

                            int ln =  gptr->data[j][0];
                            int rn =  gptr->data[j][1];
                            if ((ln == p_edge->node[1]) || (rn == p_edge->node[1]))
                            {

                                if (local_neighbor_counter == p_est->neighbor_id)
                                    //if (reservoir_sampling(p_edge->neighbor_counter) == 1)
                                {

                                    p_second_edge->id = j;
                                    p_second_edge->node[0] =  (p_edge->node[1] == ln) ? ln : rn;
                                    p_second_edge->node[1] =  (p_edge->node[1] == ln) ? rn : ln;
                                    p_second_edge->p = p_edge->neighbor_counter;

                                    p_est->status = STATUS_CLOSE;
                                    break;
                                }
                                local_neighbor_counter ++;
                            }
#if 0
                            else if ((ln == p_edge->node[0]) || (rn == p_edge->node[0]))
                            {

                                if (local_neighbor_counter == p_est->neighbor_id)
                                    //if (reservoir_sampling(p_edge->neighbor_counter) == 1)
                                {

                                    p_second_edge->id = j;
                                    p_second_edge->node[0] =  (p_edge->node[0] == ln) ? ln : rn;
                                    p_second_edge->node[1] =  (p_edge->node[0] == ln) ? rn : ln;
                                    p_second_edge->p = p_edge->neighbor_counter;

                                    p_est->status = STATUS_CLOSE;
                                    break;
                                }
                                local_neighbor_counter ++;
                            }
#endif
                        }
                        break;
                    }
                    case STATUS_CLOSE:
                    {


                        p_est->expecation = 0;
#if 1
                        int node = p_edge->node[p_est->neighbor_flag ] ;
                        //printf("node %d -> %d\n", node, p_second_edge->node[1]);

                        for (int j = csr->rpao[node]; j < csr->rpao[node + 1]; j++)
                        {
                            //printf("%d\n", csr->ciao[j]);
                            if (csr->ciao[j] == p_second_edge->node[1])
                            {
                                success_counter ++;
                                p_est->expecation = p_edge->p *  p_second_edge->p;
                                break;
                            }
                        }
                       
#else
                        for (int j = p_second_edge->id + 1; j < edgeNum; j++)
                        {
                            int ln =  gptr->data[j][0];
                            int rn =  gptr->data[j][1];
                            if (((ln == p_edge->node[p_est->neighbor_flag]) && (rn == p_second_edge->node[1])) ||
                                ((rn == p_edge->node[p_est->neighbor_flag]) && (ln == p_second_edge->node[1])))

                            {
                                success_counter ++;
                                p_est->expecation = p_edge->p *  p_second_edge->p;
                                break;
                            }
                        }
#endif
                        p_est->status = STATUS_END;
                        break;

                    }
                    }

                }
            }
            double result = 0;
            unsigned int  failed_counter = 0;
            for (int i = 0; i < est_num; i++)
            {

                estimator * p_est    = &local_estimator[i];
                sample_edge * p_edge    = &p_est->first_edge;
                sample_edge * p_second_edge = &p_est->second_edge;
#if 0
                DEBUG_PRINTF("id %d nc %d %d status %d [%d %d] [%d %d]\n",
                             p_edge->id,
                             p_edge->neighbor_counter,
                             p_est->neighbor_id,
                             p_second_edge->id,
                             p_edge->node[0],
                             p_edge->node[1],
                             p_second_edge->node[0],
                             p_second_edge->node[1]);
#endif
                result     += p_est->expecation;
            }
            DEBUG_PRINTF("result %lf @ %d with %d,total %d \n", (double(result) / est_num),
                         repeat, est_num,
                         success_counter);
        }
    }
//    for (unsigned i = 0; i < 10; ++i)
//        std::cout << ui(mt) << std::endl;

    free(gptr);
    free(csr);

    return 0;
}

