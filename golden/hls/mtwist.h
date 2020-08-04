// Mersenne Twister adapted for HLS synthesis
// author Wojtek Fedorko
// loosely based on/started from o C++ implementation from
// http://www.bedaux.net/mtrand
// and the original C implementation from M. Matsumoto and T. Nishimura (http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html)


#ifndef MTWIST_H
#define MTWIST_H

#include "ap_int.h"
#include "hls_stream.h"
#include <stdint.h>

const ap_uint<1> ap_uint_1_1=1;

const unsigned int ui_n=624;
const unsigned int ui_cycles=8;
const unsigned int ui_cycle_len=ui_n/ui_cycles;



const ap_uint<10> n=ui_n;
const ap_uint<10> cycles=ui_cycles;
const ap_uint<10> cycle_len=ui_cycle_len;

const ap_uint<10> m=397;

const ap_uint<5> MT_u = 11;
const ap_uint<5> MT_s = 7;
const ap_uint<5> MT_t = 15;
const ap_uint<5> MT_l = 18;


const ap_uint<32> MT_a = 0x9908B0DFUL;
const ap_uint<32> MT_b = 0x9D2C5680UL;
const ap_uint<32> MT_c = 0xEFC60000UL;

const ap_uint<32> MT_f = 1812433253UL;

const ap_uint<32> upper_mask=0x80000000UL;
const ap_uint<32> lower_mask=0x7FFFFFFFUL;


typedef hls::stream<ap_uint<32> > hlss;

//ap_uint<32> state[ui_n];



inline ap_uint<32> twiddle(ap_uint<32> u, ap_uint<32> v) {

  if(v.test(0))
    return (((u & upper_mask) | (v & lower_mask)) >> 1) ^ MT_a;
  return (((u & upper_mask) | (v & lower_mask)) >> 1);
//WTF: test lower bit using the ap_uint member function
}

void seed(ap_uint<32> s);

//main core
void mtwist_core(bool init, ap_uint<32> seed, ap_uint<32> stream_length, hlss& out_stream); 

#endif
