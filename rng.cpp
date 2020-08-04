
#include "approximation.h"

int rng_test(CSR* csr, int num)
{
    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0))  + counter);

    counter ++;
    int edgeNum   = csr ->edgeNum;
    boost::random::uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);
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
