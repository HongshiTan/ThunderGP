
#include "approximation.h"



extern prng mt;

int approximation_motifs_scheme_2(estimator *p_est, Graph* gptr, int edgeNum , int id)
{
    static int counter = 0;
    counter ++;
    memset(p_est, 0, sizeof(estimator));
    boost::random::uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);

    sample_edge * p_first    = &p_est->first_edge;
    //sample_edge * p_second    = &p_est->second_edge;


    int start_edge = static_cast<int>(ui(mt));
    p_est->first_edge.node[0] = gptr->data[start_edge][0];
    p_est->first_edge.node[1] = gptr->data[start_edge][1];
    p_est->first_edge.id = start_edge;
    p_first->update_counter ++;
    p_est->neighbor_counter = 0;
    p_est->status = 1;
    int temp_neighbor_counter = 0;



    for (int i = start_edge + 1; i < edgeNum; i++)
    {
        int adjacent_flag = 0;
        int ln = gptr->data[i][0];
        int rn = gptr->data[i][1];

#if 1
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
           temp_neighbor_counter ++;
        }
    }
    p_est->status = 3;
    if (p_est->status == 3)
    {
        p_est->neighbor_counter = temp_neighbor_counter;
        p_est->expecation = (double)(temp_neighbor_counter);
    }
    else
    {
        p_est->expecation = 0;
    }
    return  0;
}


int approximation_motifs_scheme_1(estimator *p_est, Graph* gptr, int edgeNum, int id)
{
    prng lmt;
    static int counter = 0;
    
    unsigned long seed = id + static_cast<unsigned int>(std::time(0)) + (unsigned long)(p_est) + counter;
    //DEBUG_PRINTF("seed %d \n",seed);
    lmt.seed(seed);
    counter ++;
    memset(p_est, 0, sizeof(estimator));
  

    sample_edge * p_first    = &p_est->first_edge;
    //sample_edge * p_second    = &p_est->second_edge;

    p_est->status = 0;
    p_est->expecation = 0;

    int start_index = 0;//edgeNum / 200 * (id % 100) + id ;
    int temp_neighbor_counter  = 0;
    for (int i = start_index; i < edgeNum; i++)
    {
        if (reservoir_sampling(i + 1, lmt) )
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
#if 0
            if ((ln == p_first->node[0]))
            {
                p_est->neighbor_flag = 1;
                adjacent_flag = 1;
            }
#endif
#if 1
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
                temp_neighbor_counter ++;
#if 0
                if (reservoir_sampling(temp_neighbor_counter, mt))
                {
                    p_est->second_edge.id = i;
                    p_est->neighbor_counter = temp_neighbor_counter;
                    p_est->status = 3;
                }
                else
                {
                    p_est->neighbor_counter ++;
                }
#endif
            }
        }
    }
    p_est->status = 3;
    if (p_est->status == 3)
    {
        p_est->expecation = ((double)temp_neighbor_counter);
    }
    else
    {
        p_est->expecation = 0;
    }
    return  0;
}