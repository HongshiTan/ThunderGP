
#include "approximation.h"


extern prng mt;

int approximation_triangle_scheme_2(estimator *p_est, Graph* gptr, int edgeNum, int id)
{
   
    static int counter = 0;
    
    counter ++;
    memset(p_est, 0, sizeof(estimator));
    boost::random::uniform_int_distribution<mpz_int> first_sample(0, edgeNum - 1);

    sample_edge * p_first    = &p_est->first_edge;
    //sample_edge * p_second    = &p_est->second_edge;


    int start_edge = static_cast<int>(first_sample(mt));
    p_est->first_edge.node[0] = gptr->data[start_edge][0];
    p_est->first_edge.node[1] = gptr->data[start_edge][1];
    p_est->first_edge.id = start_edge;
    p_first->update_counter ++;
    p_est->neighbor_counter = 0;
    p_est->status = 1;
    int tmp_neighbor_counter = 0;



    for (int i = start_edge + 1; i < edgeNum; i++)
    {
        int adjacent_flag = 0;
        int ln = gptr->data[i][0];
        int rn = gptr->data[i][1];

#if 1
        if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
        {
            //p_est->neighbor_flag = 1;
            adjacent_flag = 1;
        }
#endif
#if 1
        if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
        {
            //p_est->neighbor_flag = 0;
            adjacent_flag = 1;
        }
#endif
        if (adjacent_flag == 1)
        {
            tmp_neighbor_counter ++;
        }
    }
    p_est->tmp_neighbor_counter = tmp_neighbor_counter;
    if (p_est->tmp_neighbor_counter == 0)
    {
        return 0;
    }
    uniform_int_distribution<mpz_int> second_sample(0, p_est->tmp_neighbor_counter  - 1);
    int l2_edge = static_cast<int>(second_sample(mt));
    if (l2_edge != 0)
    {
        //    DEBUG_PRINTF("%d \n",l2_edge);
    }
    p_est->tmp_rng = l2_edge;
    //return 0;

    tmp_neighbor_counter = 0;
    for (int i = start_edge + 1; i < edgeNum; i++)
    {
        int adjacent_flag = 0;
        int ln = gptr->data[i][0];
        int rn = gptr->data[i][1];

#if 1
        if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
        {
            //p_est->neighbor_flag = 1;
            adjacent_flag = 1;
        }
#endif
#if 1
        if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
        {
            //p_est->neighbor_flag = 0;
            adjacent_flag = 1;
        }
#endif
        if (adjacent_flag == 1)
        {
          
            if (tmp_neighbor_counter == l2_edge)
            {
                p_est->second_edge.node[0] = gptr->data[i][0];
                p_est->second_edge.node[1] = gptr->data[i][1];
                p_est->second_edge.id = i;
                p_est->status = 2;
                break;
            }
            tmp_neighbor_counter ++;
        }
    }
    if (p_est->status != 2)
    {
        p_est->expecation = 0;
        return 0;
    }
    int node1, node2;
    if (p_est->second_edge.node[0] == p_est->first_edge.node[0])
    {
        node1 = p_est->second_edge.node[1];
        node2 = p_est->first_edge.node[1];
    }
    else if (p_est->second_edge.node[0] == p_est->first_edge.node[1])
    {
        node1 = p_est->second_edge.node[1];
        node2 = p_est->first_edge.node[0];
    }
    else if (p_est->second_edge.node[1] == p_est->first_edge.node[0])
    {
        node1 = p_est->second_edge.node[0];
        node2 = p_est->first_edge.node[1];
    }
    else if (p_est->second_edge.node[1] == p_est->first_edge.node[1])
    {
        node1 = p_est->second_edge.node[0];
        node2 = p_est->first_edge.node[0];
    }
    else
    {
        DEBUG_PRINTF("l2 %d %d \n", p_est->tmp_neighbor_counter,
                     p_est->tmp_rng);
        DEBUG_PRINTF("error connect : (%d %d) (%d %d)\n",
                     p_est->first_edge.node[0],
                     p_est->first_edge.node[1],
                     p_est->second_edge.node[0],
                     p_est->second_edge.node[1]);
    }
    p_est->expecation = 0;
    for (int i = p_est->second_edge.id + 1; i < edgeNum; i++)
    {
        int ln = gptr->data[i][0];
        int rn = gptr->data[i][1];
        int match_flag = 0;
        if ((ln == node1) && (rn == node2))
        {
            match_flag = 1;
        }
        else if ((ln == node2) && (rn == node1))
        {
            match_flag = 1;
        }
        if (match_flag == 1)
        {
            p_est->status = 3;
            p_est->expecation = (double)(p_est->tmp_neighbor_counter - p_est->tmp_rng + 1);
            return 1;
        }
    }
    return  0;
}

int approximation_triangle_scheme_1(estimator *p_est, Graph* gptr, int edgeNum, int id)
{
    prng lmt;
    prng lmt2;
    static int counter = 0;
    unsigned long seed = id + static_cast<unsigned int>(std::time(0)) + (unsigned long)(p_est) + counter;
    //DEBUG_PRINTF("seed %d \n",seed);
    lmt.seed(seed);
    lmt2.seed(seed + 10);
    counter ++;
    memset(p_est, 0, sizeof(estimator));

    sample_edge * p_first    = &p_est->first_edge;
    //sample_edge * p_second    = &p_est->second_edge;

    p_est->status = 0;
    p_est->expecation = 0;

    int start_index = 0;//edgeNum / 200 * (id % 100) + id ;
    int temp_neighbor_counter  = 0;
    int node1, node2;
    for (int i = start_index; i < edgeNum; i++)
    {
        if (reservoir_sampling(i + 1, lmt) )
        {
            p_est->first_edge.node[0] = gptr->data[i][0];
            p_est->first_edge.node[1] = gptr->data[i][1];
            p_first->update_counter ++;
            temp_neighbor_counter = 0;

            p_est->status = 1;
            p_est->first_edge.id = i;
        }
        else
        {
            int ln = gptr->data[i][0];
            int rn = gptr->data[i][1];


            int adjacent_flag = 0;

#if 1
            if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
            {
                adjacent_flag = 1;
            }
#endif
#if 1
            if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
            {
                adjacent_flag = 1;
            }
#endif
            if (adjacent_flag == 1)
            {
                temp_neighbor_counter ++;
#if 1
                if (reservoir_sampling(temp_neighbor_counter, lmt2))
                {
                    p_est->second_edge.node[0] = gptr->data[i][0];
                    p_est->second_edge.node[1] = gptr->data[i][1];
                    p_est->second_edge.id = i;
                    p_est->neighbor_counter = 0;

                    if (p_est->second_edge.node[0] == p_est->first_edge.node[0])
                    {
                        node1 = p_est->second_edge.node[1];
                        node2 = p_est->first_edge.node[1];
                    }
                    else if (p_est->second_edge.node[0] == p_est->first_edge.node[1])
                    {
                        node1 = p_est->second_edge.node[1];
                        node2 = p_est->first_edge.node[0];
                    }
                    else if (p_est->second_edge.node[1] == p_est->first_edge.node[0])
                    {
                        node1 = p_est->second_edge.node[0];
                        node2 = p_est->first_edge.node[1];
                    }
                    else if (p_est->second_edge.node[1] == p_est->first_edge.node[1])
                    {
                        node1 = p_est->second_edge.node[0];
                        node2 = p_est->first_edge.node[0];
                    }
                    else
                    {
                        DEBUG_PRINTF("l2 %d %d \n", p_est->tmp_neighbor_counter,
                                     p_est->tmp_rng);
                        DEBUG_PRINTF("error connect : (%d %d) (%d %d)\n",
                                     p_est->first_edge.node[0],
                                     p_est->first_edge.node[1],
                                     p_est->second_edge.node[0],
                                     p_est->second_edge.node[1]);
                    }
                    p_est->status = 2;
                }
                else
                {
                    if (p_est->status == 2)
                    {

                        //p_est->neighbor_counter ++;
                        int match_flag = 0;
                        if ((ln == node1) && (rn == node2))
                        {
                            match_flag = 1;
                        }
                        else if ((ln == node2) && (rn == node1))
                        {
                            match_flag = 1;
                        }
                        if (match_flag == 1)
                        {
                            p_est->neighbor_counter = temp_neighbor_counter + 1;
                            p_est->status = 3;

                            //break;
                        }
                    }
                }
#endif
            }
        }
    }
    if (p_est->status == 3)
    {
        p_est->expecation = ((double)p_est->neighbor_counter);
    }
    else
    {
        p_est->expecation = 0;
    }
    return  0;
}

