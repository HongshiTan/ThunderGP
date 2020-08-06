#include <stdio.h>
#include <string.h>
#include <ap_int.h>
#include "mtwist.h"

ap_uint<32> state[ui_cycles][ui_cycle_len];


ap_uint<10> p;




void seed(ap_uint<32> s) {  // init by 32 bit seed

#pragma HLS ARRAY_PARTITION variable=state complete dim=1


    state[0][0] = s;// & 0xFFFFFFFFUL; // for > 32 bit machines
    //state[0] = s;
    for (ap_uint<10> i = 1; i < n; ++i) {
        ap_uint<10> i_subarr = i % cycles;
        ap_uint<10> i_index = i / cycles;

        ap_uint<10> i_m1_subarr = (i - 1) % cycles;
        ap_uint<10> i_m1_index = (i - 1) / cycles;


        state[i_subarr][i_index] = MT_f * (state[i_m1_subarr][i_m1_index] ^ (state[i_m1_subarr][i_m1_index] >> 30)) + i;
        //state[i] = MT_f * (state[i-1] ^ (state[i-1] >> 30)) + i;
        // see Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier
        // in the previous versions, MSBs of the seed affect only MSBs of the array state
        // 2002/01/09 modified by Makoto Matsumoto
        //state[i] &= 0xFFFFFFFFUL; // for > 32 bit machines
    }
    //p = n; // force gen_state() to be called for next random number
}

void mtwist_core (bool init, ap_uint<32> seed_val, ap_uint<32> stream_length, hlss& out_stream)//, uint16_t strm_len)
{

    if (init) {
        seed(seed_val);
        p = 0;
    }

    ap_uint<32> x;
    ap_uint<32> temp_state;
    ap_uint<32> temp_state_next;

    ap_uint<10> q = p;
    ap_uint<10> q_subarr = q % cycles;
    ap_uint<10> q_index = q / cycles;

    bool p_wrap_next = (p == n - 2);
    ap_uint<10> p_next = (p + 1);
    if (p == (n - 1))
        p_next = 0;
    ap_uint<10> p_subarr = p % cycles;
    ap_uint<10> p_index = p / cycles;


    ap_uint<10> index_cacheA;
    ap_uint<10> index_cacheA_subarr;
    ap_uint<10> index_cacheA_index;



    ap_uint<10> index_cacheC;
    ap_uint<10> index_cacheC_subarr;
    ap_uint<10> index_cacheC_index;

    ap_uint<32> state_cacheA;//=state[m];
    ap_uint<32> state_cacheB;//=state[0];
    ap_uint<32> state_cacheC;//=state[1];

    state_cacheB = state[p_subarr][p_index];

    //state_cacheB=state[p];

    if (p < (n - m)) {
        index_cacheA = p + m;
        index_cacheC = p + 1;
    } else if (p < (n - 1) ) {
        index_cacheA = p + m - n;
        index_cacheC = p + 1;
    } else {
        index_cacheA = m - 1;
        index_cacheC = 0;
        //p=0;
    }
    index_cacheA_subarr = index_cacheA % cycles;
    index_cacheA_index = index_cacheA / cycles;
    index_cacheC_subarr = index_cacheC % cycles;
    index_cacheC_index = index_cacheC / cycles;

    state_cacheA = state[index_cacheA_subarr][index_cacheA_index];
    state_cacheC = state[index_cacheC_subarr][index_cacheC_index];

    temp_state_next = state_cacheA ^ twiddle(state_cacheB, state_cacheC);

    ap_uint<10> index_cacheA_next = index_cacheA + 1;
    if (index_cacheA == (n - 1))
        index_cacheA_next = 0;
    ap_uint<10> index_cacheC_next = index_cacheC + 1;
    if (index_cacheC == (n - 1))
        index_cacheC_next = 0;


    bool index_cacheA_wrap_next = index_cacheA_next == n - 1;
    bool index_cacheC_wrap_next = index_cacheC_next == n - 1;


    bool stop_next = false;
    if (stream_length == 0)
        stop_next = true;

    ap_uint<32> rand_index = 0;

rand_compute: while (true) { //for(ap_uint<32> rand_index=0;rand_index<stream_length; rand_index++){
#pragma HLS DEPENDENCE variable=state array inter false

#pragma HLS PIPELINE

        if (stop_next) {
            break;
        }
        stop_next = (rand_index == (stream_length - 1));

        rand_index++;

        temp_state = temp_state_next;
        x = state_cacheB;//state[p];

        q = p;
        p = p_next;
        p_subarr = p % cycles;
        p_index = p / cycles;


        state_cacheB = state_cacheC;

        index_cacheA = index_cacheA_next;
        index_cacheC = index_cacheC_next;



        index_cacheA_subarr = index_cacheA % cycles;
        index_cacheA_index = index_cacheA / cycles;

        index_cacheC_index = index_cacheC / cycles;


        //access: 0 4 2
        //access: 1 5 3
        //access: 2 6 4
        //access: 3 7 5
        //access: 4 0 6
        //access: 5 1 7
        //access: 6 2 0
        //access: 7 3 1

        switch (index_cacheA_subarr) {
        case 0:
            state_cacheA = state[0][index_cacheA_index];
            state_cacheC = state[4][index_cacheC_index];
            state[2][q_index] = temp_state;
            break;
        case 1:
            state_cacheA = state[1][index_cacheA_index];
            state_cacheC = state[5][index_cacheC_index];
            state[3][q_index] = temp_state;
            break;
        case 2:
            state_cacheA = state[2][index_cacheA_index];
            state_cacheC = state[6][index_cacheC_index];
            state[4][q_index] = temp_state;
            break;
        case 3:
            state_cacheA = state[3][index_cacheA_index];
            state_cacheC = state[7][index_cacheC_index];
            state[5][q_index] = temp_state;
            break;
        case 4:
            state_cacheA = state[4][index_cacheA_index];
            state_cacheC = state[0][index_cacheC_index];
            state[6][q_index] = temp_state;
            break;
        case 5:
            state_cacheA = state[5][index_cacheA_index];
            state_cacheC = state[1][index_cacheC_index];
            state[7][q_index] = temp_state;
            break;
        case 6:
            state_cacheA = state[6][index_cacheA_index];
            state_cacheC = state[2][index_cacheC_index];
            state[0][q_index] = temp_state;
            break;
        case 7:
            state_cacheA = state[7][index_cacheA_index];
            state_cacheC = state[3][index_cacheC_index];
            state[1][q_index] = temp_state;
            break;
        }

        temp_state_next = state_cacheA ^ twiddle(state_cacheB, state_cacheC);


        q_subarr = p_subarr;
        q_index = p_index;

        if (p_wrap_next) {
            p_next = 0;
        } else {
            p_next = p + 1;
        }

        p_wrap_next = (p == n - 2);

        if (index_cacheA_wrap_next) {
            index_cacheA_next = 0;
        } else {
            index_cacheA_next++;
        }
        index_cacheA_wrap_next = (index_cacheA == n - 2);


        if (index_cacheC_wrap_next) {
            index_cacheC_next = 0;
        } else {
            index_cacheC_next++;//=index_cacheC_next_next;//++;//=index_cacheC_incr2;
        }
        index_cacheC_wrap_next = (index_cacheC == n - 2);

        //if(p==n-1){
        //  p_next=0;
        //} else {
        //  p_next=p+1;
        //}

        x ^= (x >> MT_u);
        x ^= (x << MT_s) & MT_b;
        x ^= (x << MT_t) & MT_c;
        out_stream.write(x ^ (x >> MT_l));
    }
}



#define DATA_SIZE 2048
//define internal buffer max size
#define BURSTBUFFERSIZE 256

//TRIPCOUNT identifiers
const unsigned int c_size_min = 1;
const unsigned int c_size_max = BURSTBUFFERSIZE;
const unsigned int c_chunk_sz = DATA_SIZE;

extern "C" {
    void hls_rng(int *a, int size, int inc_value) {


        // Map pointer a to AXI4-master interface for global memory access
#pragma HLS INTERFACE m_axi port=a  offset=slave bundle=gmem max_read_burst_length=256 max_write_burst_length=256
        // We also need to map a and return to a bundled axilite slave interface
#pragma HLS INTERFACE s_axilite port=a bundle=control
#pragma HLS INTERFACE s_axilite port=size bundle=control
#pragma HLS INTERFACE s_axilite port=inc_value bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

        hls::stream< ap_uint<32> >  rn_stream;
        int burstbuffer[BURSTBUFFERSIZE];
        ap_uint<32> seed = 0x123;


#pragma HLS DATAFLOW

        mtwist_core(true, seed, size, rn_stream);

        //Per iteration of this loop perform BURSTBUFFERSIZE vector addition
        for (int i = 0; i < size; i += BURSTBUFFERSIZE) {
#pragma HLS LOOP_TRIPCOUNT min=c_size_min*c_size_min max=c_chunk_sz*c_chunk_sz/(c_size_max*c_size_max)
            int chunk_size = BURSTBUFFERSIZE;
            //boundary checks
            if ((i + BURSTBUFFERSIZE) > size)
                chunk_size = size - i;


loop1:      for (int j = 0; j < chunk_size; j++) {
#pragma HLS LOOP_TRIPCOUNT min=c_size_max max=c_chunk_sz
#pragma HLS PIPELINE II=1
                ap_uint<32> value = rn_stream.read();
                burstbuffer[j] = value;
            }
loop2:      for (int j = 0; j < chunk_size; j++) {
#pragma HLS LOOP_TRIPCOUNT min=c_size_max max=c_chunk_sz
#pragma HLS PIPELINE II=1
                a[i + j] = burstbuffer[j];
            }
        }
    }
}
