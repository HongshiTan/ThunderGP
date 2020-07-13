
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

#define  MAX_NUM_ESTIMATOR  (1000 * 1000 * 10)

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
                           100000, 200000, 500000, 1000 * 1000,
                           1000 * 1000 * 2, 1000 * 1000 * 5
                         };
    for (int k = 0; k < ARRAY_SIZE(est_num_map); k ++)
    {
        est_num = est_num_map[k];
        for (int repeat = 0; repeat < 1; repeat ++)
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
            int node_2_count = 0;
            for (int i = 0; i < est_num; i++)
            {
                estimator  * p_est =  &local_estimator[i];
                sample_edge * p_edge = &p_est->second_edge;

                int first_node = local_estimator[i].first_edge.node[1];
                int first_degree = csr->rpao[first_node + 1] - csr->rpao[first_node];


                int second_node = local_estimator[i].first_edge.node[0];
                int second_degree = csr->rpao[second_node + 1] - csr->rpao[second_node];
                second_degree = 0;
                while (1)
                {
                    int rn = static_cast<int>(ui(mt));
                    int id = int((double (rn) * ratio) * (first_degree + second_degree));
                    //printf("id %d rn %d degree %d %d \n",id,rn, first_degree,second_degree);

                    if (id > first_degree)
                    {
                        node_2_count ++;
                        int sub_id = id - first_degree;
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
                        int sub_id = id;
                        p_edge->id = sub_id;
                        p_edge->node[0] = first_node;
                        p_edge->node[1] = csr->ciao[csr->rpao[first_node] + sub_id];
                        p_edge->p = first_degree + second_degree;
                        break;
                    }

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
                    if (csr->ciao[j] == second_node && second_node != first_node)
                    {
                        p_est->expecation =  p_est->second_edge.p;
                        //printf("[%d-%d-%d], %d\n",first_node,p_est->first_edge.node[1], second_node,p_est->expecation);
                        result += p_est->expecation;

                        break;
                    }
                }
                if (p_est->expecation == 0)
                {
                    failed_counter ++;
                }
            }
            DEBUG_PRINTF("result %lf @ %d with %d,total %d - node_2_count %d \n", (double(result) / est_num),
                         repeat, est_num,
                         est_num - failed_counter, node_2_count);
        }
    }
//    for (unsigned i = 0; i < 10; ++i)
//        std::cout << ui(mt) << std::endl;

    free(gptr);
    free(csr);

    return 0;
}

