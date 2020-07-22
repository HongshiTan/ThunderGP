
#include "approximation.h"


typedef struct
{
    unsigned int first_edge_id;
    unsigned int src_node;
    unsigned int dst_node;
    unsigned int src_deg;
    unsigned int dst_deg;
    int c;
    unsigned int second_edge_rng;
    unsigned int tmp_deg;
    unsigned int tmp_deg_flag;


    unsigned int second_edge_id;
    unsigned int src_node2;
    unsigned int dst_node2;
    unsigned int node1;
    unsigned int node2;
} simple_est;

simple_est  localest[MAX_NUM_ESTIMATOR];

int triangle_count(Graph* gptr, CSR* csr, int num)
{
    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0))  + counter);

    counter ++;
    int edgeNum   = csr ->edgeNum;
    uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);
    for (int i = 0; i < num; i ++)
    {
        memset(&localest[i], 0, sizeof(simple_est));
        localest[i].first_edge_id = static_cast<int>(ui(mt));
        localest[i].src_node = gptr->data[localest[i].first_edge_id ][0];
        localest[i].dst_node = gptr->data[localest[i].first_edge_id ][1];
    }


    for (int i = 0; i < edgeNum; i++)
    {

        unsigned int src = gptr->data[i][0];
        unsigned int dst = gptr->data[i][1];
        for (int k = 0; k < num; k ++)
        {
            if (i < localest[k].first_edge_id)
            {
                if ((src == localest[k].src_node) || ((dst == localest[k].src_node)))
                {
                    localest[k].src_deg ++;
                }
                if ((src == localest[k].dst_node) || (dst == localest[k].dst_node))
                {
                    localest[k].dst_deg ++;
                }
            }
        }
    }
    for (int k = 0; k < num; k++)
    {
        int src_deg = csr->rpao[localest[k].src_node + 1] - csr->rpao[localest[k].src_node];
        int dst_deg = csr->rpao[localest[k].dst_node + 1] - csr->rpao[localest[k].dst_node];
        localest[k].c = src_deg + dst_deg  - localest[k].src_deg - localest[k].dst_deg;

        //DEBUG_PRINTF("%d %d %d %d \n", src_deg , dst_deg , localest[k].src_deg , localest[k].dst_deg );

        if (localest[k].c > 0)
        {
            uniform_int_distribution<mpz_int> ui(1, localest[k].c);
            //DEBUG_PRINTF("c %d \n", localest[k].c );
            localest[k].second_edge_rng = static_cast<int>(ui(mt));
            if (localest[k].second_edge_rng <= (src_deg -  localest[k].src_deg) )
            {
                int idx = csr->rpao[localest[k].src_node] + localest[k].second_edge_rng;
                //DEBUG_PRINTF("idx %d \n", idx );
                localest[k].dst_node2 = csr->ciao[idx];
                localest[k].src_node2 = localest[k].src_node;
                localest[k].tmp_deg_flag = 0;
                localest[k].tmp_deg = idx;
            }
            else
            {
                int idx = csr->rpao[localest[k].src_node] + localest[k].second_edge_rng
                          - (src_deg -  localest[k].src_deg);

                localest[k].dst_node2 = csr->ciao[idx];
                localest[k].src_node2 = localest[k].dst_node;
                localest[k].tmp_deg_flag = 1;
                localest[k].tmp_deg = idx;
            }
        }
        else if (localest[k].c < 0)
        {
            localest[k].tmp_deg_flag = 2;
            DEBUG_PRINTF("error %d @%d \n", localest[k].c, k);
        }
        localest[k].src_deg = 0;
        localest[k].dst_deg = 0;
    }

    for (int i = 0; i < edgeNum; i++)
    {

        unsigned int src = gptr->data[i][0];
        unsigned int dst = gptr->data[i][1];
        for (int k = 0; k < num; k ++)
        {
            if ((src == localest[k].src_node) || ((dst == localest[k].src_node)))
            {
                localest[k].src_deg ++;
            }
            if ((src == localest[k].dst_node) || (dst == localest[k].dst_node))
            {
                localest[k].dst_deg ++;
            }
            if (localest[k].tmp_deg_flag == 0)
            {
                if (localest[k].src_deg == localest[k].tmp_deg)
                {
                    localest[k].second_edge_id = k;
                }
            }
            else if (localest[k].tmp_deg_flag == 1)
            {
                if (localest[k].dst_deg == localest[k].tmp_deg)
                {
                    localest[k].second_edge_id = k;
                }
            }
        }
    }
    for (int k = 0; k < num; k++)
    {
        if (localest[k].tmp_deg_flag == 0)
        {
            localest[k].node1 = localest[k].dst_node;
            localest[k].node2 = localest[k].dst_node2;
        }
        else if (localest[k].tmp_deg_flag == 1)
        {
            localest[k].node1 = localest[k].src_node;
            localest[k].node2 = localest[k].dst_node2;
        }
    }

    double result = 0;
    for (int i = 0; i < edgeNum; i++)
    {
        unsigned int ln = gptr->data[i][0];
        unsigned int rn = gptr->data[i][1];
        for (int k = 0; k < num; k++)
        {
            if ((i > localest[k].second_edge_id) && (localest[k].tmp_deg_flag < 2))
            {
#if 1
                if ((ln == localest[k].node1) && (rn == localest[k].node2))
                {
                    result += localest[k].c;

                }
                else if ((rn == localest[k].node1) && (ln == localest[k].node2))
                {
                    result += localest[k].c;
                }
#endif
            }
        }
    }
    DEBUG_PRINTF("[result] %d %f\n", num, (result * edgeNum / (double)(num)));
    return 0;
}

int approximation_triangle_scheme_2(estimator *p_est, Graph* gptr, CSR* csr, int id)
{
    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0)) + (unsigned long)(p_est) + counter);

    counter ++;
    memset(p_est, 0, sizeof(estimator));
    int edgeNum   = csr ->edgeNum;
    uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);

    sample_edge * p_first    = &p_est->first_edge;
    sample_edge * p_second   = &p_est->second_edge;


    int start_edge = static_cast<int>(ui(mt));
    p_est->first_edge.node[0] = gptr->data[start_edge][0];
    p_est->first_edge.node[1] = gptr->data[start_edge][1];
    p_est->first_edge.id = start_edge;
    p_first->update_counter ++;
    p_est->neighbor_counter = 0;
    p_est->status = 1;
    int temp_neighbor_counter = 0;


    return 0;
}


int approximation_triangle_scheme_1(estimator *p_est, Graph* gptr, CSR* csr, int id)
{

    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0)) + (unsigned long)(p_est) + counter);
    counter ++;
    memset(p_est, 0, sizeof(estimator));
    int edgeNum   = csr ->edgeNum;

    sample_edge * p_first    = &p_est->first_edge;
    sample_edge * p_second    = &p_est->second_edge;

    p_est->status = 0;
    p_est->expecation = 0;

    int start_index = 0;//edgeNum / 200 * (id % 100) + id ;
    int temp_neighbor_counter  = 0;
    for (int i = start_index; i < edgeNum; i++)
    {
        if (reservoir_sampling(i + 1, mt) )
        {
            p_est->first_edge.node[0] = gptr->data[i][0];
            p_est->first_edge.node[1] = gptr->data[i][1];
            p_first->update_counter ++;
            temp_neighbor_counter = 0;
            p_est->neighbor_counter = 0;
            p_est->status = 1;
            p_est->first_edge.id = i;
        }
        else
        {
            int ln = gptr->data[i][0];
            int rn = gptr->data[i][1];


            int adjacent_flag = 0;
#if 1
            if ((ln == p_first->node[0]))
            {
                p_est->neighbor_flag = 1;
                adjacent_flag = 1;
            }
#endif
#if 0
            if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
            {
                p_est->neighbor_flag = 1;
                adjacent_flag = 1;
            }
#endif
#if 0
            if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
            {
                p_est->neighbor_flag = 0;
                adjacent_flag = 1;
            }
#endif
            if (adjacent_flag == 1)
            {
                temp_neighbor_counter ++;
#if 1
                if (reservoir_sampling(temp_neighbor_counter, mt))
                {
                    p_est->second_edge.id = i;
                    p_est->neighbor_counter = temp_neighbor_counter;
                    p_est->status = 2;
                    p_second->update_counter ++;

                    if ( p_est->neighbor_flag == 0)
                    {
                        p_second->node[0] =  (p_first->node[1] == ln) ? ln : rn;
                        p_second->node[1] =  (p_first->node[1] == ln) ? rn : ln;
                    }
                    else
                    {
                        p_second->node[0] =  (p_first->node[0] == ln) ? ln : rn;
                        p_second->node[1] =  (p_first->node[0] == ln) ? rn : ln;
                    }
                }
                else
                {
                    p_est->neighbor_counter ++;
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
#endif
            }
        }
    }
    if (p_est->status == 3)
    {
        p_est->expecation = temp_neighbor_counter  * edgeNum;
    }
    else
    {
        p_est->expecation = 0;
    }
    return  0;
}
