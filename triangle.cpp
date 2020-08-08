
#include "approximation.h"

#include <list>

#define S3_DEBUG(fmt,...)    ;


extern prng mt;

typedef struct
{
    int id;
    int deg;
    int v;
} l2_sample;

typedef struct
{
    int id;
    int r2_id;
    int node[2];
} check_sample;




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


        if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
        {
            adjacent_flag = 1;
        }
        if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
        {
            adjacent_flag = 1;
        }

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
    p_est->neighbor_counter = p_est->tmp_neighbor_counter;
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


        if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
        {
            //p_est->neighbor_flag = 1;
            adjacent_flag = 1;
        }
        if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
        {
            //p_est->neighbor_flag = 0;
            adjacent_flag = 1;
        }

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
            p_est->success = 1;
            p_est->expecation = (double)(p_est->tmp_neighbor_counter - p_est->tmp_rng + 1);
            return 1;
        }
    }
    p_est->success = 0;
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

            if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
            {
                adjacent_flag = 1;
            }
            if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
            {
                adjacent_flag = 1;
            }
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


int approximation_triangle_scheme_4(estimator *p_est, Graph* gptr, int edgeNum, int id)
{
    pcg32_random_t lmt;
    pcg32_random_t inner;


    static int counter = 0;
    //unsigned long seed = id + static_cast<unsigned int>(std::time(0)) + (unsigned long)(p_est) + counter;
    pcg32_srandom_r(&lmt, 0x853c49e6748fea9bULL,  0xda3e39cb94b95bdbULL + counter);
    pcg32_srandom_r(&inner, 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL + counter + id);
    //DEBUG_PRINTF("seed %d \n",seed);
    //lmt.seed(seed);
    //lmt2.seed(seed + 10);

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
        if (pcg_reservoir_sampling(i + 1, &lmt) )
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

            if ((ln == p_first->node[0]) || (rn == p_first->node[0]))
            {
                adjacent_flag ++;
            }
            if ((ln == p_first->node[1]) || (rn == p_first->node[1]))
            {
                adjacent_flag ++;
            }
            if (adjacent_flag != 0 )
            {
                if(adjacent_flag > 1)
                {
                    DEBUG_PRINTF("bug\n");
                }
                temp_neighbor_counter ++;
#if 1
                if (pcg_reservoir_sampling(temp_neighbor_counter, &inner))
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
        p_est->success = 1;
        p_est->expecation = ((double)p_est->neighbor_counter);
    }
    else
    {
        p_est->success = 0;
        p_est->expecation = 0;
    }
    return  0;
}


int approximation_triangle_scheme_3(estimator *g_est, Graph* gptr, int edgeNum , int est_id)
{
    estimator sub[SUB_EST];
    int id[SUB_EST];
    static int counter = 0;
    counter ++;
    int * vertex_deg = (int *)memalign(64 * 1024, gptr->vertexNum * 4);

    int * tmp_deg = (int *)memalign(64 * 1024, gptr->vertexNum * 4);
    memset(tmp_deg, 0, gptr->vertexNum * 4);
    std::list<l2_sample> tmp_deg_list[SUB_EST];


    int * close_check = (int *)memalign(64 * 1024, gptr->vertexNum * 4);
    memset(close_check, 0, gptr->vertexNum * 4);
    std::list<check_sample> close_check_list[SUB_EST * 2];



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
    memset(vertex_deg, 0, gptr->vertexNum * 4);
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
    int current_list_count = 0;
    pthread_mutex_lock(&lock);
    for (int i = 0; i < SUB_EST; i++)
    {
        int ln = sub[i].first_edge.node[0];
        int rn = sub[i].first_edge.node[1];
        int a = vertex_deg[ln] - sub[i].first_edge.temp_deg[0];
        int b = vertex_deg[rn] - sub[i].first_edge.temp_deg[1];
        sub[i].neighbor_counter = a + b;
        //g_est->expecation += sub[i].neighbor_counter;
        if (sub[i].neighbor_counter > 0)
        {
            uniform_int_distribution<mpz_int> second_sample(1, sub[i].neighbor_counter);
            
            int rng = static_cast<int>(second_sample(mt));
            sub[i].tmp_rng = rng;
            int selected_deg;
            int selected_v;

            if (rng > a)
            {
                //rn
                selected_deg = rng - a  + sub[i].first_edge.temp_deg[1];
                selected_v = rn;
            }
            else
            {
                //ln
                selected_deg = rng + sub[i].first_edge.temp_deg[0];
                selected_v = ln;
            }
            l2_sample r2;
            r2.id  = i;
            r2.deg = selected_deg;
            r2.v = selected_v;

            if (tmp_deg[selected_v] == 0)
            {
                tmp_deg[selected_v] = current_list_count;
                current_list_count ++;
            }
            else
            {
                S3_DEBUG("Add to tail of list %d with size %ld v %d \n",
                         tmp_deg[selected_v],
                         tmp_deg_list[tmp_deg[selected_v]].size(), selected_v);
                S3_DEBUG("##################################################################\n")
            }
            S3_DEBUG("info: a %d b %d rng %d,[%d-%d] %d, c %d\n",
                     a, b, rng, ln, rn, selected_v, current_list_count);
            S3_DEBUG("deg:[%d %d] tmp:[%d %d] - %d\n",
                     vertex_deg[ln],
                     vertex_deg[rn],
                     sub[i].first_edge.temp_deg[0],
                     sub[i].first_edge.temp_deg[1],
                     selected_deg);

            tmp_deg_list[tmp_deg[selected_v]].push_back(r2);

            //g_est->success ++;
        }
    }
    pthread_mutex_unlock(&lock);

    memset(vertex_deg, 0, gptr->vertexNum * 4);
    for (int i = 0; i < edgeNum; i++)
    {
        int ln = gptr->data[i][0];
        int rn = gptr->data[i][1];
        vertex_deg[ln] ++;
        vertex_deg[rn] ++;
        if (tmp_deg[ln] != 0)
        {
            for (auto item : tmp_deg_list[tmp_deg[ln]])
            {
                int s_v = item.v;
                int id = item.id;
                if (vertex_deg[s_v] == item.deg)
                {
                    sub[id].second_edge.node[0] = ln;
                    sub[id].second_edge.node[1] = rn;
                    sub[id].second_edge.id = i;
                    sub[id].status = 2;
                    S3_DEBUG("[%d] r1 %d r2 %d \n", id, sub[id].first_edge.id, i);
                }
            }
        }
        if (tmp_deg[rn] != 0)
        {
            for (auto item : tmp_deg_list[tmp_deg[rn]])
            {
                int s_v = item.v;
                int id = item.id;
                if (vertex_deg[s_v] == item.deg)
                {
                    sub[id].second_edge.node[0] = ln;
                    sub[id].second_edge.node[1] = rn;
                    sub[id].second_edge.id = i;
                    sub[id].status = 2;
                    S3_DEBUG("[%d] r1 %d r2 %d \n", id, sub[id].first_edge.id, i);
                }
            }
        }
    }

    // close check
    int current_close_list_count = 0;
    for (int i = 0; i < SUB_EST; i++)
    {
        int node1, node2;
        estimator *p_est = &sub[i];
        if (p_est->status == 2)
        {

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
                continue;
            }

            //g_est->expecation += sub[i].neighbor_counter;
            check_sample close;

            close.id = i;
            close.r2_id = p_est->second_edge.id;
            close.node[0] = node1;
            close.node[1] = node2;


            if (close_check[node1] == 0)
            {
                close_check[node1] = current_close_list_count;
                current_close_list_count ++;
            }
            close_check_list[close_check[node1]].push_back(close);

            if (close_check[node2] == 0)
            {
                close_check[node2] = current_close_list_count;
                current_close_list_count ++;
            }
            close_check_list[close_check[node2]].push_back(close);
        }
    }
    double result = 0;
    int success = 0;
    for (int i = 0; i < edgeNum; i++)
    {
        int ln = gptr->data[i][0];
        int rn = gptr->data[i][1];
        if (close_check[ln] != 0)
        {
            for (auto item : close_check_list[close_check[ln]])
            {
                int k = item.r2_id;
                int id = item.id;
                int node1 = item.node[0];
                int node2 = item.node[1];
                if ( i > k)
                {
                    int sampled_flag = 0;
                    if ((node1 == ln) && (node2 == rn))
                    {
                        sampled_flag = 1;
                    }
                    else if ((node2 == ln) && (node1 == rn))
                    {
                        sampled_flag = 1;
                    }
                    if (sampled_flag == 1)
                    {
                        if (sub[id].status == 2)
                        {
                            result += sub[id].neighbor_counter;
                            sub[id].status = 3;
                            success ++;
                        }
                    }
                }
            }
        }
        if (close_check[rn] != 0)
        {
            for (auto item : close_check_list[close_check[rn]])
            {
                int k = item.r2_id;
                int id = item.id;
                int node1 = item.node[0];
                int node2 = item.node[1];
                if ( i > k)
                {
                    int sampled_flag = 0;
                    if ((node1 == ln) && (node2 == rn))
                    {
                        sampled_flag = 1;
                    }
                    else if ((node2 == ln) && (node1 == rn))
                    {
                        sampled_flag = 1;
                    }
                    if (sampled_flag == 1)
                    {
                        if (sub[id].status == 2)
                        {
                            result += sub[id].neighbor_counter;
                            sub[id].status = 3;
                            success ++;
                        }
                    }
                }
            }
        }
    }


    g_est->status = 3;
    g_est->expecation = (double)result;
    g_est->success = success;
    free(vertex_deg);
    free(tmp_deg);
    free(close_check);
    S3_DEBUG("------------------------------------------------------------\n")
    return  0;
}