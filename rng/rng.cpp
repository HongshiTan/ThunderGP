
#include <stdint.h>
#include "approximation.h"








uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_rng_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_rng_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_rng_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}


uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound)
{
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.

    uint32_rng_t threshold = -bound % bound;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In 
    // practice, bounds are typically small and only a tiny amount of the range
    // is eliminated.
    for (;;) {
        uint32_rng_t r = pcg32_random_r(rng);
        if (r >= threshold)
            return r % bound;
    }
}


int reservoir_sampling(int n, prng &lmt)
{
    boost::random::bernoulli_distribution<double> p(1.0 / (n));
    pthread_mutex_lock(&lock);
    bool res = p(mt);
    pthread_mutex_unlock(&lock);
    //DEBUG_PRINTF("n %d, res %d\n", n, res);
    if (res)
    {
        return 1;
    }
    return 0;
}




int pcg_reservoir_sampling_2(int n, pcg32_random_t* rng)
{
    uint32_t value = pcg32_random_r(rng);
    double p = pow(2,32) / ((double)n); 
    //uint32_t value = pcg32_boundedrand_r(rng, n<<4);
    //boost::random::bernoulli_distribution<double> p(1.0 / (n));
    //pthread_mutex_lock(&lock);
    //bool res = p(mt);
    //pthread_mutex_unlock(&lock);
    //DEBUG_PRINTF("n %d, res %d\n", n, res);
    if (value < p )
    {
        return 1;
    }
    return 0;
}


int pcg_reservoir_sampling(int n, pcg32_random_t* rng)
{
    uint32_t value = pcg32_random_r(rng);
    double p = (double (pow(2,32))) / ((double)n); 
    //uint32_t value = pcg32_boundedrand_r(rng, n << 4);
    //boost::random::bernoulli_distribution<double> p(1.0 / (n));
    //pthread_mutex_lock(&lock);
    //bool res = p(mt);
    //pthread_mutex_unlock(&lock);
    //DEBUG_PRINTF("n %d, res %d\n", n, res);
    //if (value < (1 << 4))
    if (value < p )
    {
        return 1;
    }
    return 0;
}

//2987624

int rng_test(int edgeNum, int num)
{
    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0))  + counter);

    counter ++;
    boost::random::uniform_int_distribution<mpz_int> ui(0, edgeNum - 1);
    for (int i = 0; i < num; i ++)
    {
        int loc = static_cast<int>(ui(mt));
        DEBUG_PRINTF("%d %d \n", loc , edgeNum);
    }
    return 0;
}

int rng_res_test(int edgeNum, int num)
{
    prng mt;
    static int counter = 0;
    mt.seed(static_cast<unsigned int>(std::time(0))  + counter);

    counter ++;
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

        DEBUG_PRINTF("%d %d \n", loc , edgeNum);
    }
    return 0;
}

int rng_pcg_res_test(int edgeNum, int num)
{
    pcg32_random_t mt;
    static int counter = 0;
    pcg32_srandom_r(&mt, 0x853c49e6748fea9bULL,  0xda3e39cb94b95bdbULL );

    counter ++;
    for (int i = 0; i < num; i ++)
    {
        int loc = 0;
        for (int k = 0; k < edgeNum; k ++)
        {
            if (pcg_reservoir_sampling(k + 1, &mt) )
            {
                loc = k;
            }
        }

        DEBUG_PRINTF("%d %d \n", loc , edgeNum);
    }
    return 0;
}
