
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
    unsigned int node[2];
    unsigned int p;
} sample_edge;


typedef struct {
    unsigned int expecation;
    sample_edge first_edge;
    sample_edge second_edge;
} estimator;

#define  MAX_NUM_ESTIMATOR  (1000 * 1000)

estimator local_estimator[MAX_NUM_ESTIMATOR];

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

    std::cout << std::hex << std::showbase;
    std::cout << (mpz_int(1) << RNG_RANGE_MAX) << std::endl;

    Graph* gptr = createGraph(gName, mode);
    CSR* csr    = new CSR(*gptr);

    int vertexNum = csr ->vertexNum;
    int edgeNum   = csr ->edgeNum;


//ARRAY_SIZE(local_estimator)
    int est_num_map [] = { 10, 100, 1000, 2000,
                           5000, 10000, 20000, 50000,
                           100000, 200000, 500000
                         };
    for (int k = 0; k < 1; k ++)
    {
        est_num = est_num_map[k];
        for (int repeat = 0; repeat < 10; repeat ++)
        {
            double ratio = 1.0 / (double(1 << RNG_RANGE_MAX));
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
            }
            for (int i = 0; i < est_num; i++)
            {
                estimator  * p_est =  &local_estimator[i];
                sample_edge * p_edge = &p_est->second_edge;
                int rn = static_cast<int>(ui(mt));
                int first_node = local_estimator[i].first_edge.node[1];
                int first_degree = csr->rpao[first_node + 1] - csr->rpao[first_node];


                int second_node = local_estimator[i].first_edge.node[0];
                int second_degree = csr->rpao[second_node + 1] - csr->rpao[second_node];

                int id = int((double (rn)) * ratio * (first_degree + second_degree));
                if (id >= first_node)
                {
                    int sub_id = id - first_node;
                    p_edge->id = sub_id;
                    p_edge->node[0] = second_node;
                    p_edge->node[1] = csr->ciao[csr->rpao[second_node] + sub_id];
                    p_edge->p = first_degree + second_degree;
                    int tmp = p_est->first_edge.node[0];
                    p_est->first_edge.node[0] = p_est->first_edge.node[1];
                    p_est->first_edge.node[1] = tmp;
                }
                else
                {
                    int sub_id = first_node;
                    p_edge->id = sub_id;
                    p_edge->node[0] = first_node;
                    p_edge->node[1] = csr->ciao[csr->rpao[first_node] + sub_id];
                    p_edge->p = first_degree + second_degree;
                }

            }
            unsigned long long result = 0;
            unsigned int  failed_counter = 0;
            for (int i = 0; i < est_num; i++)
            {
                estimator  * p_est =  &local_estimator[i];
                int first_node = p_est->first_edge.node[0];
                int second_node = p_est->second_edge.node[1];

                p_est->expecation = 0;
                for (int j = csr->rpao[first_node]; j < csr->rpao[first_node + 1]; j ++ )
                {
                    if (csr->ciao[j] == second_node)
                    {
                        p_est->expecation = p_est->first_edge.p * p_est->second_edge.p;
                        result += p_est->expecation;

                        break;
                    }
                }
                if (p_est->expecation == 0)
                {
                    failed_counter ++;
                }
            }
            DEBUG_PRINTF("result %llu @ %d with %d,failed_counter %d\n", result / est_num, repeat, est_num, failed_counter);
        }
    }
//    for (unsigned i = 0; i < 10; ++i)
//        std::cout << ui(mt) << std::endl;

    free(gptr);
    free(csr);

    return 0;
}

