/* Effect: change sample rate  Copyright (c) 2008,12 robs@users.sourceforge.net
 *
 * Copyright (C) 2012 lvqcl - modifications
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#define onehalf 0.5f

static __inline __m128 _mm_add_horz_ss(__m128 x)
{
    __m128 y;
    y = _mm_shuffle_ps(x, x, _MM_SHUFFLE(2, 3, 0, 1));
    x = _mm_add_ps(x, y);
    y = _mm_movehl_ps(y, x);
    x = _mm_add_ss(x, y);
    return x;
}

#ifdef _MSC_VER
#pragma warning(disable : 4305)
#endif

static _MM_ALIGN16 const sox_sample_t half_fir_coefs_8[] = {
  0.3115465451887802, -0.08734497241282892, 0.03681452335604365,
  -0.01518925831569441, 0.005454118437408876, -0.001564400922162005,
  0.0003181701445034203, -3.48001341225749e-5,
};

static _MM_ALIGN16 const sox_sample_t half_fir_coefs_9[] = {
  0.3122703613711853, -0.08922155288172305, 0.03913974805854332,
  -0.01725059723447163, 0.006858970092378141, -0.002304518467568703,
  0.0006096426006051062, -0.0001132393923815236, 1.119795386287666e-5,
};

static _MM_ALIGN16 const sox_sample_t half_fir_coefs_10[] = {
  0.3128545521327376, -0.09075671986104322, 0.04109637155154835,
  -0.01906629512749895, 0.008184039342054333, -0.0030766775017262,
  0.0009639607022414314, -0.0002358552746579827, 4.025184282444155e-5,
  -3.629779111541012e-6,
};

static _MM_ALIGN16 const sox_sample_t half_fir_coefs_11[] = {
  0.3133358837508807, -0.09203588680609488, 0.04276515428384758,
  -0.02067356614745591, 0.00942253142371517, -0.003856330993895144,
  0.001363470684892284, -0.0003987400965541919, 9.058629923971627e-5,
  -1.428553070915318e-5, 1.183455238783835e-6,
};

static _MM_ALIGN16 const sox_sample_t half_fir_coefs_12[] = {
  0.3137392991811407, -0.0931182192961332, 0.0442050575271454,
  -0.02210391200618091, 0.01057473015666001, -0.00462766983973885,
  0.001793630226239453, -0.0005961819959665878, 0.0001631475979359577,
  -3.45557865639653e-5, 5.06188341942088e-6, -3.877010943315563e-7,
};

static _MM_ALIGN16 const sox_sample_t half_fir_coefs_13[] = {
  0.3140822554324578, -0.0940458550886253, 0.04545990399121566,
  -0.02338339450796002, 0.01164429409071052, -0.005380686021429845,
  0.002242915773871009, -0.000822047600000082, 0.0002572510962395222,
  -6.607320708956279e-5, 1.309926399120154e-5, -1.790719575255006e-6,
  1.27504961098836e-7,
};

#ifdef _MSC_VER
#pragma warning(default : 4305)
#endif


#define HSUM0(cf) \
    sum = _mm_set_ss(0.5f); \
    x1 = _mm_loadu_ps(input  ); \
    x2 = _mm_loadu_ps(input+4); \
    x3 = _mm_loadu_ps(input-4); \
    x4 = _mm_loadu_ps(input-8); \
    x5 = _mm_load_ps(half_fir_coefs_##cf); \
    sum = _mm_mul_ps(sum, x1); \
    x1 = _mm_shuffle_ps(x1, x2, _MM_SHUFFLE(3, 1, 3, 1)); \
    x3 = _mm_shuffle_ps(x3, x4, _MM_SHUFFLE(1, 3, 1, 3)); \
    x1 = _mm_add_ps(x1, x3); \
    x1 = _mm_mul_ps(x1, x5); \
    sum = _mm_add_ps(sum, x1)

#define HSUM(cf, n) \
    x1 = _mm_loadu_ps(input+2*n  ); \
    x2 = _mm_loadu_ps(input+2*n+4); \
    x3 = _mm_loadu_ps(input-2*n-4); \
    x4 = _mm_loadu_ps(input-2*n-8); \
    x5 = _mm_load_ps(half_fir_coefs_##cf + n); \
    x1 = _mm_shuffle_ps(x1, x2, _MM_SHUFFLE(3, 1, 3, 1)); \
    x3 = _mm_shuffle_ps(x3, x4, _MM_SHUFFLE(1, 3, 1, 3)); \
    x1 = _mm_add_ps(x1, x3); \
    x1 = _mm_mul_ps(x1, x5); \
    sum = _mm_add_ps(sum, x1)

#define HSUM1(cf, n) \
    x1 = _mm_loadu_ps(input+2*n); \
    x3 = _mm_loadu_ps(input-2*n-4); \
    x1 = _mm_shuffle_ps(x1, x3, _MM_SHUFFLE(1, 3, 3, 1)); \
    x5 = _mm_loadl_pi(x5, (__m64*)(half_fir_coefs_##cf + n)); \
    x5 = _mm_movelh_ps(x5, x5); \
    x1 = _mm_mul_ps(x1, x5); \
    sum = _mm_add_ps(sum, x1)

#define HSUM1A(cf, n) \
    x1 = _mm_load_ss(input+2*n+1); \
    x3 = _mm_load_ss(input-2*n-1); \
    x1 = _mm_add_ss(x1, x3); \
    x5 = _mm_load_ss(half_fir_coefs_##cf + n); \
    x1 = _mm_mul_ss(x1, x5); \
    sum = _mm_add_ss(sum, x1)


/* Down-sample by a factor of 2 using a FIR with odd length (LEN).*/
/* Input must be preceded and followed by LEN >> 1 samples. */

static void h8(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre; /* = stage_read_p(p) */
  int i, num_out = (stage_occupancy(p) + 1) / 2;
  sox_sample_t * const output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2)
  {
    __m128 sum, x1, x2, x3, x4, x5;
    //sox_sample_t sum = input[0] * onehalf;
    //sum += (input[-1 ] + input[1 ]) * half_fir_coefs_8[0];
    //sum += (input[-3 ] + input[3 ]) * half_fir_coefs_8[1];
    //sum += (input[-5 ] + input[5 ]) * half_fir_coefs_8[2];
    //sum += (input[-7 ] + input[7 ]) * half_fir_coefs_8[3];
    HSUM0(8);
    //sum += (input[-9 ] + input[9 ]) * half_fir_coefs_8[4];
    //sum += (input[-11] + input[11]) * half_fir_coefs_8[5];
    //sum += (input[-13] + input[13]) * half_fir_coefs_8[6];
    //sum += (input[-15] + input[15]) * half_fir_coefs_8[7];
    HSUM(8, 4);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, 2 * num_out, NULL);
}

static void h10(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_out = (stage_occupancy(p) + 1) / 2;
  sox_sample_t * const output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2)
  {
    __m128 sum, x1, x2, x3, x4, x5;

    HSUM0(10);
    HSUM(10, 4);

    //sum += (input[-17] + input[17]) * half_fir_coefs_10[8];
    //sum += (input[-19] + input[19]) * half_fir_coefs_10[9];
    HSUM1(10, 8);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, 2 * num_out, NULL);
}

static void h12(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_out = (stage_occupancy(p) + 1) / 2;
  sox_sample_t * const output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2)
  {
    __m128 sum, x1, x2, x3, x4, x5;

    HSUM0(12);
    HSUM(12, 4);
    HSUM(12, 8);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, 2 * num_out, NULL);
}

static void h9(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_out = (stage_occupancy(p) + 1) / 2;
  sox_sample_t * const output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2)
  {
    __m128 sum, x1, x2, x3, x4, x5;

    HSUM0(9);
    HSUM(9, 4);
    //sum += (input[-17] + input[17]) * half_fir_coefs_9[8];
    HSUM1A(9, 8);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, 2 * num_out, NULL);
}

static void h11(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_out = (stage_occupancy(p) + 1) / 2;
  sox_sample_t * const output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2)
  {
    __m128 sum, x1, x2, x3, x4, x5;

    HSUM0(11);
    HSUM(11, 4);
    HSUM1(11, 8);
    HSUM1A(11, 10);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, 2 * num_out, NULL);
}

static void h13(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_out = (stage_occupancy(p) + 1) / 2;
  sox_sample_t * const output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2)
  {
    __m128 sum, x1, x2, x3, x4, x5;

    HSUM0(12);
    HSUM(12, 4);
    HSUM(12, 8);
    HSUM1A(12, 12);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, 2 * num_out, NULL);
}

#ifdef _MSC_VER
#pragma warning(disable : 4305)
#endif

static struct {int num_coefs; stage_fn_t fn; float att;} const half_firs[] = {
  { 8, h8 , 136.51},
  { 9, h9 , 152.32},
  {10, h10, 168.07},
  {11, h11, 183.78},
  {12, h12, 199.44},
  {13, h13, 212.75},
};

#ifdef _MSC_VER
#pragma warning(default : 4305)
#endif

#define SUMM4(n) \
      xmm1 = _mm_loadu_ps(at+n); \
      xmm1 = _mm_mul_ps(xmm1, _mm_load_ps(coefs+n)); \
      sum = _mm_add_ps(sum, xmm1)

#define SUMM4p1z(j, k) \
        xmm3 = _mm_loadu_ps(in+j); \
        xmm2 = _mm_load_ps(coefs+k); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+4)); \
        xmm2 = _mm_mul_ps(xmm2, xmm3); \
        sum = _mm_add_ps(sum, xmm2)

#define SUMM4p2z(j, k) \
        xmm3 = _mm_loadu_ps(in+j); \
        xmm2 = _mm_load_ps(coefs+k); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+4)); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+8)); \
        xmm2 = _mm_mul_ps(xmm2, xmm3); \
        sum = _mm_add_ps(sum, xmm2)

#define SUMM4p3z(j, k) \
        xmm3 = _mm_loadu_ps(in+j); \
        xmm2 = _mm_load_ps(coefs+k); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+4)); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+8)); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+12)); \
        xmm2 = _mm_mul_ps(xmm2, xmm3); \
        sum = _mm_add_ps(sum, xmm2)

#define SUMM8(n) \
      xmm1 = _mm_loadu_ps(at+n); \
      xmm2 = _mm_loadu_ps(at+n+4); \
      xmm1 = _mm_mul_ps(xmm1, _mm_load_ps(coefs+n)); \
      xmm2 = _mm_mul_ps(xmm2, _mm_load_ps(coefs+n+4)); \
      sum = _mm_add_ps(sum, xmm1); \
      sum = _mm_add_ps(sum, xmm2)

#define SUMM8p1z(j, k) \
        xmm3 = _mm_loadu_ps(in+j); \
        xmm5 = _mm_loadu_ps(in+j+4); \
        xmm2 = _mm_load_ps(coefs+k); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+4)); \
        xmm4 = _mm_load_ps(coefs+k+8); \
        xmm4 = _mm_mul_ps(xmm4, xx1); \
        xmm4 = _mm_add_ps(xmm4, _mm_load_ps(coefs+k+12)); \
        xmm2 = _mm_mul_ps(xmm2, xmm3); \
        sum = _mm_add_ps(sum, xmm2); \
        xmm4 = _mm_mul_ps(xmm4, xmm5); \
        sum = _mm_add_ps(sum, xmm4)

#define SUMM8p2z(j, k) \
        xmm3 = _mm_loadu_ps(in+j); \
        xmm5 = _mm_loadu_ps(in+j+4); \
        xmm2 = _mm_load_ps(coefs+k); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+4)); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+8)); \
        xmm4 = _mm_load_ps(coefs+k+12); \
        xmm4 = _mm_mul_ps(xmm4, xx1); \
        xmm4 = _mm_add_ps(xmm4, _mm_load_ps(coefs+k+16)); \
        xmm4 = _mm_mul_ps(xmm4, xx1); \
        xmm4 = _mm_add_ps(xmm4, _mm_load_ps(coefs+k+20)); \
        xmm2 = _mm_mul_ps(xmm2, xmm3); \
        sum = _mm_add_ps(sum, xmm2); \
        xmm4 = _mm_mul_ps(xmm4, xmm5); \
        sum = _mm_add_ps(sum, xmm4)

#define SUMM8p3z(j, k) \
        xmm3 = _mm_loadu_ps(in+j); \
        xmm5 = _mm_loadu_ps(in+j+4); \
        xmm2 = _mm_load_ps(coefs+k); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+4)); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+8)); \
        xmm2 = _mm_mul_ps(xmm2, xx1); \
        xmm2 = _mm_add_ps(xmm2, _mm_load_ps(coefs+k+12)); \
        xmm4 = _mm_load_ps(coefs+k+16); \
        xmm4 = _mm_mul_ps(xmm4, xx1); \
        xmm4 = _mm_add_ps(xmm4, _mm_load_ps(coefs+k+20)); \
        xmm4 = _mm_mul_ps(xmm4, xx1); \
        xmm4 = _mm_add_ps(xmm4, _mm_load_ps(coefs+k+24)); \
        xmm4 = _mm_mul_ps(xmm4, xx1); \
        xmm4 = _mm_add_ps(xmm4, _mm_load_ps(coefs+k+28)); \
        xmm2 = _mm_mul_ps(xmm2, xmm3); \
        sum = _mm_add_ps(sum, xmm2); \
        xmm4 = _mm_mul_ps(xmm4, xmm5); \
        sum = _mm_add_ps(sum, xmm4)


#define SUMM4p1(j) SUMM4p1z(j, (j)*2)
#define SUMM4p2(j) SUMM4p2z(j, (j)*3)
#define SUMM4p3(j) SUMM4p3z(j, (j)*4)

#define SUMM8p1(j) SUMM8p1z(j, (j)*2)
#define SUMM8p2(j) SUMM8p2z(j, (j)*3)
#define SUMM8p3(j) SUMM8p3z(j, (j)*4)

#define SWITCHER(func4, func8) \
      case  0: break; \
      case  4: func4(0); break; \
      case  8: func8(0); break; \
      case 12: func8(0); func4(8); break; \
      case 16: func8(0); func8(8); break; \
      case 20: func8(0); func8(8); func4(16); break; \
      case 24: func8(0); func8(8); func8(16); break; \
      case 28: func8(0); func8(8); func8(16); func4(24); break; \
      default: func8(0); func8(8); func8(16); func8(24); j = 32;


/* Resample using a non-interpolated poly-phase FIR with length LEN.*/
/* Input must be followed by LEN-1 samples. */

static void vpoly0(stage_t * p, fifo_t * output_fifo)
{
  int num_in = stage_occupancy(p);
  if (num_in > 0)
  {
    int i;
    const int max_num_out = (num_in * p->L - p->at.parts.integer + p->step.parts.integer - 1)/p->step.parts.integer;
    sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
    sox_sample_t * const output = (sox_sample_t *)fifo_reserve(output_fifo, max_num_out);
    sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;
    div_t divided2;

    assert(p->n%4 == 0);
    for (i = 0; p->at.parts.integer < num_in * p->L; ++i, p->at.parts.integer += p->step.parts.integer) {
      div_t divided = div1(p->at.parts.integer, p->L);
      sox_sample_t const * const at = input + divided.quot;
      sox_sample_t const * const coefs = coef0 + p->n*divided.rem;
      int j = 0;
      __m128 sum, xmm1, xmm2;
      sum = _mm_setzero_ps();
      switch(p->n) {
        SWITCHER(SUMM4, SUMM8);
        while (j < p->n) {
          SUMM4(j);
          j+=4;
        }
      }
      sum = _mm_add_horz_ss(sum);
      _mm_store_ss(output+i, sum);
    }
    assert(max_num_out - i >= 0);
    //fifo_trim_by(output_fifo, max_num_out - i);
    divided2 = div1(p->at.parts.integer, p->L);
    fifo_read(&p->fifo, divided2.quot, NULL);
    p->at.parts.integer = divided2.rem;
  }
}


/* Resample using an interpolated poly-phase FIR with length LEN.*/
/* Input must be followed by LEN-1 samples. */

static void vpoly1(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;

  assert(p->n%4 == 0);
  if (p->use_hi_prec_clock) {
    hi_prec_clock_t at = p->at.hi_prec_clock;
    for (i = 0; (int)at < num_in; ++i, at += p->step.hi_prec_clock) {
      sox_sample_t const * const in = input + (int)at;
      hi_prec_clock_t fraction = at - (int)at;
      int phase = fraction * (1 << p->phase_bits);
      sox_sample_t const * const coefs = coef0 + p->n*phase*2;

      sox_sample_t x = fraction * (1 << p->phase_bits) - phase;
      int j = 0;
      int k = 0;
      __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
      sum = _mm_setzero_ps();
      xx1 = _mm_set1_ps(x);
      switch(p->n) {
        SWITCHER(SUMM4p1, SUMM8p1);
        while (j < p->n) {
          SUMM4p1z(j,k);
          j+=4; k+=8;
        }
      }
      sum = _mm_add_horz_ss(sum);
      _mm_store_ss(output+i, sum);
    }
    fifo_read(&p->fifo, (int)at, NULL);
    p->at.hi_prec_clock = at - (int)at;
  } else
  {
    for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
      sox_sample_t const * const in = input + p->at.parts.integer;
      uint32_t fraction = p->at.parts.fraction;
      int phase = fraction >> (32 - p->phase_bits); /* high-order bits */
      sox_sample_t const * const coefs = coef0 + p->n*phase*2;

      sox_sample_t x = (sox_sample_t) (fraction << p->phase_bits) * (1 / MULT32);
      int j = 0;
      int k = 0;
      __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
      sum = _mm_setzero_ps();
      xx1 = _mm_set1_ps(x);
      switch(p->n) {
        SWITCHER(SUMM4p1, SUMM8p1);
        while (j < p->n) {
          SUMM4p1z(j,k);
          j+=4; k+=8;
        }
      }
      sum = _mm_add_horz_ss(sum);
      _mm_store_ss(output+i, sum);
    }
    fifo_read(&p->fifo, p->at.parts.integer, NULL);
    p->at.parts.integer = 0;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
}

static void vpoly2(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;

  assert(p->n%4 == 0);
  if (p->use_hi_prec_clock) {
    hi_prec_clock_t at = p->at.hi_prec_clock;
    for (i = 0; (int)at < num_in; ++i, at += p->step.hi_prec_clock) {
      sox_sample_t const * const in = input + (int)at;
      hi_prec_clock_t fraction = at - (int)at;
      int phase = fraction * (1 << p->phase_bits);
      sox_sample_t const * const coefs = coef0 + p->n*phase*3;

      sox_sample_t x = fraction * (1 << p->phase_bits) - phase;
      int j = 0;
      int k = 0;
      __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
      sum = _mm_setzero_ps();
      xx1 = _mm_set1_ps(x);
      switch(p->n) {
        SWITCHER(SUMM4p2, SUMM8p2);
        while (j < p->n) {
          SUMM4p2z(j, k);
          j+=4; k+=12;
        }
      }
      sum = _mm_add_horz_ss(sum);
      _mm_store_ss(output+i, sum);
    }
    fifo_read(&p->fifo, (int)at, NULL);
    p->at.hi_prec_clock = at - (int)at;
  } else
  {
    for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
      sox_sample_t const * const in = input + p->at.parts.integer;
      uint32_t fraction = p->at.parts.fraction;
      int phase = fraction >> (32 - p->phase_bits); /* high-order bits */
      sox_sample_t const * const coefs = coef0 + p->n*phase*3;

      sox_sample_t x = (sox_sample_t) (fraction << p->phase_bits) * (1 / MULT32);
      int j = 0;
      int k = 0;
      __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
      sum = _mm_setzero_ps();
      xx1 = _mm_set1_ps(x);
      switch(p->n) {
        SWITCHER(SUMM4p2, SUMM8p2);
        while (j < p->n) {
          SUMM4p2z(j, k);
          j+=4; k+=12;
        }
      }
      sum = _mm_add_horz_ss(sum);
      _mm_store_ss(output+i, sum);
    }
    fifo_read(&p->fifo, p->at.parts.integer, NULL);
    p->at.parts.integer = 0;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
}

static void vpoly3(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;

  assert(p->n%4 == 0);
  if (p->use_hi_prec_clock) {
    hi_prec_clock_t at = p->at.hi_prec_clock;
    for (i = 0; (int)at < num_in; ++i, at += p->step.hi_prec_clock) {
      sox_sample_t const * const in = input + (int)at;
      hi_prec_clock_t fraction = at - (int)at;
      int phase = fraction * (1 << p->phase_bits);
      sox_sample_t const * const coefs = coef0 + p->n*phase*4;

      sox_sample_t x = fraction * (1 << p->phase_bits) - phase;
      int j = 0;
      int k = 0;
      __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
      sum = _mm_setzero_ps();
      xx1 = _mm_set1_ps(x);
      switch(p->n) {
        SWITCHER(SUMM4p3, SUMM8p3);
        while (j < p->n) {
          SUMM4p3z(j, k);
          j+=4; k+=16;
        }
      }
      sum = _mm_add_horz_ss(sum);
      _mm_store_ss(output+i, sum);
    }
    fifo_read(&p->fifo, (int)at, NULL);
    p->at.hi_prec_clock = at - (int)at;
  } else
  {
    for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
      sox_sample_t const * const in = input + p->at.parts.integer;
      uint32_t fraction = p->at.parts.fraction;
      int phase = fraction >> (32 - p->phase_bits); /* high-order bits */
      sox_sample_t const * const coefs = coef0 + p->n*phase*4;

      sox_sample_t x = (sox_sample_t) (fraction << p->phase_bits) * (1 / MULT32);
      int j = 0;
      int k = 0;
      __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
      sum = _mm_setzero_ps();
      xx1 = _mm_set1_ps(x);
      switch(p->n) {
        SWITCHER(SUMM4p3, SUMM8p3);
        while (j < p->n) {
          SUMM4p3z(j, k);
          j+=4; k+=16;
        }
      }
      sum = _mm_add_horz_ss(sum);
      _mm_store_ss(output+i, sum);
    }
    fifo_read(&p->fifo, p->at.parts.integer, NULL);
    p->at.parts.integer = 0;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
}


#define U100_l 44
#define u100_l 12
#define u100_1_b 8
#define u100_2_b 6


/* Resample using a non-interpolated poly-phase FIR with length LEN.*/
/* Input must be followed by LEN-1 samples. */

static void U100_0(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;
  div_t divided2;

  for (i = 0; p->at.parts.integer < num_in * p->L; ++i, p->at.parts.integer += p->step.parts.integer) {
    div_t divided = div1(p->at.parts.integer, p->L);
    sox_sample_t const * const at = input + divided.quot;
    sox_sample_t const * const coefs = coef0 + U100_l*divided.rem;
    __m128 sum, xmm1, xmm2;
    sum = _mm_setzero_ps();

    SUMM8(0);
    SUMM8(8);
    SUMM8(16);
    SUMM8(24);
    SUMM8(32);
    SUMM4(40);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  divided2 = div1(p->at.parts.integer, p->L);
  fifo_read(&p->fifo, divided2.quot, NULL);
  p->at.parts.integer = divided2.rem;
}

static void u100_0(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;
  div_t divided2;

  for (i = 0; p->at.parts.integer < num_in * p->L; ++i, p->at.parts.integer += p->step.parts.integer) {
    div_t divided = div1(p->at.parts.integer, p->L);
    sox_sample_t const * const at = input + divided.quot;
    sox_sample_t const * const coefs = coef0 + u100_l*divided.rem;
    __m128 sum, xmm1, xmm2;
    sum = _mm_setzero_ps();

    SUMM8(0);
    SUMM4(8);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  divided2 = div1(p->at.parts.integer, p->L);
  fifo_read(&p->fifo, divided2.quot, NULL);
  p->at.parts.integer = divided2.rem;
}


/* Resample using an interpolated poly-phase FIR with length LEN.*/
/* Input must be followed by LEN-1 samples. */

static void u100_1(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;

  for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
    sox_sample_t const * const in = input + p->at.parts.integer;
    uint32_t fraction = p->at.parts.fraction;
    int phase = fraction >> (32 - u100_1_b); /* high-order bits */
    sox_sample_t const * const coefs = coef0 + 2*u100_l*phase;

    sox_sample_t x = (sox_sample_t) (fraction << u100_1_b) * (1 / MULT32);
    __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
    sum = _mm_setzero_ps();
    xx1 = _mm_set1_ps(x);

    SUMM8p1(0);
    SUMM4p1(8);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, p->at.parts.integer, NULL);
  p->at.parts.integer = 0;

  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
}

static void u100_2(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;

  for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
    sox_sample_t const * const in = input + p->at.parts.integer;
    uint32_t fraction = p->at.parts.fraction;
    int phase = fraction >> (32 - u100_2_b); /* high-order bits */
    sox_sample_t const * const coefs = coef0 + 3*u100_l*phase;

    sox_sample_t x = (sox_sample_t) (fraction << u100_2_b) * (1 / MULT32);
    __m128 sum, xx1, xmm2, xmm3, xmm4, xmm5;
    sum = _mm_setzero_ps();
    xx1 = _mm_set1_ps(x);

    SUMM8p2(0);
    SUMM4p2(8);

    sum = _mm_add_horz_ss(sum);
    _mm_store_ss(output+i, sum);
  }
  fifo_read(&p->fifo, p->at.parts.integer, NULL);
  p->at.parts.integer = 0;

  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
}


typedef struct {float scalar; stage_fn_t fn;} poly_fir1_t;
typedef struct {float beta; poly_fir1_t interp[3];} poly_fir_t;

#ifdef _MSC_VER
#pragma warning(disable : 4305)
#endif

static poly_fir_t const poly_firs[] = {
  {-1, {{0, vpoly0}, { 7.2, vpoly1}, {5.0, vpoly2}}},
  {-1, {{0, vpoly0}, { 9.4, vpoly1}, {6.7, vpoly2}}},
  {-1, {{0, vpoly0}, {12.4, vpoly1}, {7.8, vpoly2}}},
  {-1, {{0, vpoly0}, {13.6, vpoly1}, {9.3, vpoly2}}},
  {-1, {{0, vpoly0}, {10.5, vpoly2}, {8.4, vpoly3}}},
  {-1, {{0, vpoly0}, {11.85,vpoly2}, {9.0, vpoly3}}},

  {-1, {{0, vpoly0}, { 8.0, vpoly1}, {5.3, vpoly2}}},
  {-1, {{0, vpoly0}, { 8.6, vpoly1}, {5.7, vpoly2}}},
  {-1, {{0, vpoly0}, {10.6, vpoly1}, {6.75,vpoly2}}},
  {-1, {{0, vpoly0}, {12.6, vpoly1}, {8.6, vpoly2}}},
  {-1, {{0, vpoly0}, { 9.6, vpoly2}, {7.6, vpoly3}}},
  {-1, {{0, vpoly0}, {11.4, vpoly2}, {8.65,vpoly3}}},

  {10.62, {{U100_l, U100_0}, {0, 0}, {0, 0}}},
  {11.28, {{u100_l, u100_0}, {u100_1_b, u100_1}, {u100_2_b, u100_2}}},
  {-1, {{0, vpoly0}, {   9, vpoly1}, {  6, vpoly2}}},
  {-1, {{0, vpoly0}, {  11, vpoly1}, {  7, vpoly2}}},
  {-1, {{0, vpoly0}, {  13, vpoly1}, {  8, vpoly2}}},
  {-1, {{0, vpoly0}, {  10, vpoly2}, {  8, vpoly3}}},
  {-1, {{0, vpoly0}, {  12, vpoly2}, {  9, vpoly3}}},
};

#ifdef _MSC_VER
#pragma warning(default : 4305)
#endif
