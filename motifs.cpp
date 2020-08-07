
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


        if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
        {
            p_est->neighbor_flag = 1;
            adjacent_flag = 1;
        }
        if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
        {
            p_est->neighbor_flag = 0;
            adjacent_flag = 1;
        }

        if (adjacent_flag == 1)
        {
            temp_neighbor_counter ++;
        }
    }
    p_est->status = 3;
    if (p_est->status == 3)
    {
        p_est->success  = 1;
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


            if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
            {
                p_est->neighbor_flag = 1;
                adjacent_flag = 1;
            }
            if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
            {
                p_est->neighbor_flag = 0;
                adjacent_flag = 1;
            }

            if (adjacent_flag == 1)
            {
                temp_neighbor_counter ++;
            }
        }
    }
    p_est->status = 3;
    if (p_est->status == 3)
    {
        p_est->success  = 1;
        p_est->expecation = ((double)temp_neighbor_counter);
    }
    else
    {
        p_est->expecation = 0;
    }
    return  0;
}

int approximation_motifs_scheme_3(estimator *g_est, Graph* gptr, int edgeNum , int est_id)
{
    estimator sub[SUB_EST];
    int id[SUB_EST];
    static int counter = 0;
    counter ++;
    int * vertex_deg = (int *)memalign(64 * 1024, gptr->vertexNum * 4);

    boost::random::uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);

    pthread_mutex_lock(&lock);
    for (int i = 0; i < SUB_EST; i++)
    {
        id[i] = static_cast<int>(ui(mt));
    }
    pthread_mutex_unlock(&lock);

    for (int i = 0; i < SUB_EST; i++)
    {
        for (int j = 0; j < SUB_EST - i - 1; j++)
        {
            if (id[j] > id[j + 1])
            {
                int tmpId = id[j];
                id[j] = id[j + 1];
                id[j + 1]  = tmpId;
            }
        }
    }

    for (int i = 0; i < SUB_EST; i++)
    {
        int start_edge = id[i];
        estimator *p_est = &sub[i];
        memset(p_est, 0, sizeof(estimator));

        sample_edge * p_first    = &p_est->first_edge;
        //sample_edge * p_second    = &p_est->second_edge;


        p_est->first_edge.node[0] = gptr->data[start_edge][0];
        p_est->first_edge.node[1] = gptr->data[start_edge][1];
        p_est->first_edge.id = start_edge;
        p_first->update_counter ++;
        p_est->neighbor_counter = 0;
        p_est->status = 1;
        p_est->temp_neighbor_counter = 0;
        //DEBUG_PRINTF("%d \n",id[i]);
    }
    //DEBUG_PRINTF("done \n");
    //return 0;
    int current_est_id = 0;
    for (int i = 0; i < edgeNum; i++)
    {
        int ln = gptr->data[i][0];
        int rn = gptr->data[i][1];
        vertex_deg[ln] ++;
        vertex_deg[rn] ++;
        if ((unsigned int)i == sub[current_est_id].first_edge.id)
        {
            for (int j = current_est_id ; j < SUB_EST; j ++)
            {
                if ((unsigned int)i == sub[j].first_edge.id)
                {
                    sub[j].first_edge.temp_deg[0] = vertex_deg[ln];
                    sub[j].first_edge.temp_deg[1] = vertex_deg[rn];
                }
                else
                {
                    current_est_id = j;
                    break;
                }
            }
        }

    }
    g_est->status = 3;
    g_est->expecation = 0;
    for (int i = 0; i < SUB_EST; i++)
    {
        int ln = sub[i].first_edge.node[0];
        int rn = sub[i].first_edge.node[1];
        sub[i].neighbor_counter = vertex_deg[ln] + vertex_deg[rn] - sub[i].first_edge.temp_deg[0] - sub[i].first_edge.temp_deg[1];
        g_est->expecation += sub[i].neighbor_counter;
        if (sub[i].neighbor_counter > 0)
        {
            g_est->success ++;
        }
    }
    free(vertex_deg);
    return  0;
}