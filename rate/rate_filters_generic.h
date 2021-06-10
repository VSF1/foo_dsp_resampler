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

#ifdef RATE_SINGLE_PRECISION
#define onehalf 0.5f
#else
#define onehalf 0.5
#endif


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


/* Down-sample by a factor of 2 using a FIR with odd length (LEN).*/
/* Input must be preceded and followed by LEN >> 1 samples. */

static void h8(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre; /* = stage_read_p(p) */
  int i, num_out = (stage_occupancy(p) + 1) / 2;
  sox_sample_t * const output = fifo_reserve(output_fifo, num_out);

  for (i = 0; i < num_out; ++i, input += 2)
  {
    sox_sample_t sum = input[0] * onehalf;

    sum += (input[-1 ] + input[1 ]) * half_fir_coefs_8[0];
    sum += (input[-3 ] + input[3 ]) * half_fir_coefs_8[1];
    sum += (input[-5 ] + input[5 ]) * half_fir_coefs_8[2];
    sum += (input[-7 ] + input[7 ]) * half_fir_coefs_8[3];

    sum += (input[-9 ] + input[9 ]) * half_fir_coefs_8[4];
    sum += (input[-11] + input[11]) * half_fir_coefs_8[5];
    sum += (input[-13] + input[13]) * half_fir_coefs_8[6];
    sum += (input[-15] + input[15]) * half_fir_coefs_8[7];

    output[i] = sum;
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
    sox_sample_t sum = input[0] * onehalf;

    sum += (input[-1 ] + input[1 ]) * half_fir_coefs_9[0];
    sum += (input[-3 ] + input[3 ]) * half_fir_coefs_9[1];
    sum += (input[-5 ] + input[5 ]) * half_fir_coefs_9[2];
    sum += (input[-7 ] + input[7 ]) * half_fir_coefs_9[3];

    sum += (input[-9 ] + input[9 ]) * half_fir_coefs_9[4];
    sum += (input[-11] + input[11]) * half_fir_coefs_9[5];
    sum += (input[-13] + input[13]) * half_fir_coefs_9[6];
    sum += (input[-15] + input[15]) * half_fir_coefs_9[7];

    sum += (input[-17] + input[17]) * half_fir_coefs_9[8];

    output[i] = sum;
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
    sox_sample_t sum = input[0] * onehalf;

    sum += (input[-1 ] + input[1 ]) * half_fir_coefs_10[0];
    sum += (input[-3 ] + input[3 ]) * half_fir_coefs_10[1];
    sum += (input[-5 ] + input[5 ]) * half_fir_coefs_10[2];
    sum += (input[-7 ] + input[7 ]) * half_fir_coefs_10[3];

    sum += (input[-9 ] + input[9 ]) * half_fir_coefs_10[4];
    sum += (input[-11] + input[11]) * half_fir_coefs_10[5];
    sum += (input[-13] + input[13]) * half_fir_coefs_10[6];
    sum += (input[-15] + input[15]) * half_fir_coefs_10[7];

    sum += (input[-17] + input[17]) * half_fir_coefs_10[8];
    sum += (input[-19] + input[19]) * half_fir_coefs_10[9];

    output[i] = sum;
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
    sox_sample_t sum = input[0] * onehalf;

    sum += (input[-1 ] + input[1 ]) * half_fir_coefs_11[0 ];
    sum += (input[-3 ] + input[3 ]) * half_fir_coefs_11[1 ];
    sum += (input[-5 ] + input[5 ]) * half_fir_coefs_11[2 ];
    sum += (input[-7 ] + input[7 ]) * half_fir_coefs_11[3 ];

    sum += (input[-9 ] + input[9 ]) * half_fir_coefs_11[4 ];
    sum += (input[-11] + input[11]) * half_fir_coefs_11[5 ];
    sum += (input[-13] + input[13]) * half_fir_coefs_11[6 ];
    sum += (input[-15] + input[15]) * half_fir_coefs_11[7 ];

    sum += (input[-17] + input[17]) * half_fir_coefs_11[8 ];
    sum += (input[-19] + input[19]) * half_fir_coefs_11[9 ];
    sum += (input[-21] + input[21]) * half_fir_coefs_11[10];

    output[i] = sum;
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
    sox_sample_t sum = input[0] * onehalf;

    sum += (input[-1 ] + input[1 ]) * half_fir_coefs_12[0 ];
    sum += (input[-3 ] + input[3 ]) * half_fir_coefs_12[1 ];
    sum += (input[-5 ] + input[5 ]) * half_fir_coefs_12[2 ];
    sum += (input[-7 ] + input[7 ]) * half_fir_coefs_12[3 ];

    sum += (input[-9 ] + input[9 ]) * half_fir_coefs_12[4 ];
    sum += (input[-11] + input[11]) * half_fir_coefs_12[5 ];
    sum += (input[-13] + input[13]) * half_fir_coefs_12[6 ];
    sum += (input[-15] + input[15]) * half_fir_coefs_12[7 ];

    sum += (input[-17] + input[17]) * half_fir_coefs_12[8 ];
    sum += (input[-19] + input[19]) * half_fir_coefs_12[9 ];
    sum += (input[-21] + input[21]) * half_fir_coefs_12[10];
    sum += (input[-23] + input[23]) * half_fir_coefs_12[11];

    output[i] = sum;
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
    sox_sample_t sum = input[0] * onehalf;

    sum += (input[-1 ] + input[1 ]) * half_fir_coefs_13[0 ];
    sum += (input[-3 ] + input[3 ]) * half_fir_coefs_13[1 ];
    sum += (input[-5 ] + input[5 ]) * half_fir_coefs_13[2 ];
    sum += (input[-7 ] + input[7 ]) * half_fir_coefs_13[3 ];

    sum += (input[-9 ] + input[9 ]) * half_fir_coefs_13[4 ];
    sum += (input[-11] + input[11]) * half_fir_coefs_13[5 ];
    sum += (input[-13] + input[13]) * half_fir_coefs_13[6 ];
    sum += (input[-15] + input[15]) * half_fir_coefs_13[7 ];

    sum += (input[-17] + input[17]) * half_fir_coefs_13[8 ];
    sum += (input[-19] + input[19]) * half_fir_coefs_13[9 ];
    sum += (input[-21] + input[21]) * half_fir_coefs_13[10];
    sum += (input[-23] + input[23]) * half_fir_coefs_13[11];

    sum += (input[-25] + input[25]) * half_fir_coefs_13[12];

    output[i] = sum;
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


/* Resample using a non-interpolated poly-phase FIR with length LEN.*/
/* Input must be followed by LEN-1 samples. */

static void vpoly0(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;
  div_t divided2;

  assert(p->n%4 == 0);
  for (i = 0; p->at.parts.integer < num_in * p->L; ++i, p->at.parts.integer += p->step.parts.integer) {
    div_t divided = div(p->at.parts.integer, p->L);
    sox_sample_t const * const at = input + divided.quot;
    sox_sample_t const * const coefs = coef0 + p->n*divided.rem;
    sox_sample_t sum = 0;
    int j = 0;
    int k = 0;
    while (j < p->n)
    {
      sum += coefs[k+0] * at[j+0];
      sum += coefs[k+1] * at[j+1];

      sum += coefs[k+2] * at[j+2];
      sum += coefs[k+3] * at[j+3];

      j+=4; k+=4;
    }
    output[i] = sum;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  divided2 = div(p->at.parts.integer, p->L);
  fifo_read(&p->fifo, divided2.quot, NULL);
  p->at.parts.integer = divided2.rem;
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
      sox_sample_t sum = 0;
      int j = 0;
      int k = 0;
      while (j < p->n)
      {
        sum += (coefs[k+0]*x + coefs[k+1]) * in[j+0];
        sum += (coefs[k+2]*x + coefs[k+3]) * in[j+1];

        sum += (coefs[k+4]*x + coefs[k+5]) * in[j+2];
        sum += (coefs[k+6]*x + coefs[k+7]) * in[j+3];

        j+=4; k+=8;
      }
      output[i] = sum;
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
      sox_sample_t sum = 0;
      int j = 0;
      int k = 0;
      while (j < p->n)
      {
        sum += (coefs[k+0]*x + coefs[k+1]) * in[j+0];
        sum += (coefs[k+2]*x + coefs[k+3]) * in[j+1];

        sum += (coefs[k+4]*x + coefs[k+5]) * in[j+2];
        sum += (coefs[k+6]*x + coefs[k+7]) * in[j+3];

        j+=4; k+=8;
      }
      output[i] = sum;
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
      sox_sample_t sum = 0;
      int j = 0;
      int k = 0;
      while (j < p->n)
      {
        sum += ((coefs[k+0]*x + coefs[k+1 ])*x + coefs[k+2 ]) * in[j+0];
        sum += ((coefs[k+3]*x + coefs[k+4 ])*x + coefs[k+5 ]) * in[j+1];

        sum += ((coefs[k+6]*x + coefs[k+7 ])*x + coefs[k+8 ]) * in[j+2];
        sum += ((coefs[k+9]*x + coefs[k+10])*x + coefs[k+11]) * in[j+3];

        j+=4; k+=12;
      }
      output[i] = sum;
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
      sox_sample_t sum = 0;
      int j = 0;
      int k = 0;
      while (j < p->n)
      {
        sum += ((coefs[k+0]*x + coefs[k+1 ])*x + coefs[k+2 ]) * in[j+0];
        sum += ((coefs[k+3]*x + coefs[k+4 ])*x + coefs[k+5 ]) * in[j+1];

        sum += ((coefs[k+6]*x + coefs[k+7 ])*x + coefs[k+8 ]) * in[j+2];
        sum += ((coefs[k+9]*x + coefs[k+10])*x + coefs[k+11]) * in[j+3];

        j+=4; k+=12;
      }
      output[i] = sum;
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
      sox_sample_t sum = 0;
      int j = 0;
      int k = 0;
      while (j < p->n)
      {
        sum += (((coefs[k+0 ]*x + coefs[k+1 ])*x + coefs[k+2 ])*x + coefs[k+3 ]) * in[j+0];
        sum += (((coefs[k+4 ]*x + coefs[k+5 ])*x + coefs[k+6 ])*x + coefs[k+7 ]) * in[j+1];

        sum += (((coefs[k+8 ]*x + coefs[k+9 ])*x + coefs[k+10])*x + coefs[k+11]) * in[j+2];
        sum += (((coefs[k+12]*x + coefs[k+13])*x + coefs[k+14])*x + coefs[k+15]) * in[j+3];

        j+=4; k+=16;
      }
      output[i] = sum;
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
      sox_sample_t sum = 0;
      int j = 0;
      int k = 0;
      while (j < p->n)
      {
        sum += (((coefs[k+0 ]*x + coefs[k+1 ])*x + coefs[k+2 ])*x + coefs[k+3 ]) * in[j+0];
        sum += (((coefs[k+4 ]*x + coefs[k+5 ])*x + coefs[k+6 ])*x + coefs[k+7 ]) * in[j+1];

        sum += (((coefs[k+8 ]*x + coefs[k+9 ])*x + coefs[k+10])*x + coefs[k+11]) * in[j+2];
        sum += (((coefs[k+12]*x + coefs[k+13])*x + coefs[k+14])*x + coefs[k+15]) * in[j+3];

        j+=4; k+=16;
      }
      output[i] = sum;
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
    div_t divided = div(p->at.parts.integer, p->L);
    sox_sample_t const * const at = input + divided.quot;
    sox_sample_t const * const coefs = coef0 + U100_l*divided.rem;
    sox_sample_t sum = 0;

    sum += coefs[0 ] * at[0 ];
    sum += coefs[1 ] * at[1 ];
    sum += coefs[2 ] * at[2 ];
    sum += coefs[3 ] * at[3 ];

    sum += coefs[4 ] * at[4 ];
    sum += coefs[5 ] * at[5 ];
    sum += coefs[6 ] * at[6 ];
    sum += coefs[7 ] * at[7 ];

    sum += coefs[8 ] * at[8 ];
    sum += coefs[9 ] * at[9 ];
    sum += coefs[10] * at[10];
    sum += coefs[11] * at[11];

    sum += coefs[12] * at[12];
    sum += coefs[13] * at[13];
    sum += coefs[14] * at[14];
    sum += coefs[15] * at[15];

    sum += coefs[16] * at[16];
    sum += coefs[17] * at[17];
    sum += coefs[18] * at[18];
    sum += coefs[19] * at[19];

    sum += coefs[20] * at[20];
    sum += coefs[21] * at[21];
    sum += coefs[22] * at[22];
    sum += coefs[23] * at[23];

    sum += coefs[24] * at[24];
    sum += coefs[25] * at[25];
    sum += coefs[26] * at[26];
    sum += coefs[27] * at[27];

    sum += coefs[28] * at[28];
    sum += coefs[29] * at[29];
    sum += coefs[30] * at[30];
    sum += coefs[31] * at[31];

    sum += coefs[32] * at[32];
    sum += coefs[33] * at[33];
    sum += coefs[34] * at[34];
    sum += coefs[35] * at[35];

    sum += coefs[36] * at[36];
    sum += coefs[37] * at[37];
    sum += coefs[38] * at[38];
    sum += coefs[39] * at[39];

    sum += coefs[40] * at[40];
    sum += coefs[41] * at[41];
    sum += coefs[42] * at[42];
    sum += coefs[43] * at[43];

    output[i] = sum;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  divided2 = div(p->at.parts.integer, p->L);
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
    div_t divided = div(p->at.parts.integer, p->L);
    sox_sample_t const * const at = input + divided.quot;
    sox_sample_t const * const coefs = coef0 + u100_l*divided.rem;
    sox_sample_t sum = 0;

    sum += coefs[0 ] * at[0 ];
    sum += coefs[1 ] * at[1 ];
    sum += coefs[2 ] * at[2 ];
    sum += coefs[3 ] * at[3 ];

    sum += coefs[4 ] * at[4 ];
    sum += coefs[5 ] * at[5 ];
    sum += coefs[6 ] * at[6 ];
    sum += coefs[7 ] * at[7 ];

    sum += coefs[8 ] * at[8 ];
    sum += coefs[9 ] * at[9 ];
    sum += coefs[10] * at[10];
    sum += coefs[11] * at[11];

    output[i] = sum;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  divided2 = div(p->at.parts.integer, p->L);
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
  {
    for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
      sox_sample_t const * const in = input + p->at.parts.integer;
      uint32_t fraction = p->at.parts.fraction;
      int phase = fraction >> (32 - u100_1_b); /* high-order bits */
      sox_sample_t const * const coefs = coef0 + 2*u100_l*phase;

      sox_sample_t x = (sox_sample_t) (fraction << u100_1_b) * (1 / MULT32);
      sox_sample_t sum = 0;

      sum += (coefs[0 ]*x + coefs[1 ]) * in[0 ];
      sum += (coefs[2 ]*x + coefs[3 ]) * in[1 ];
      sum += (coefs[4 ]*x + coefs[5 ]) * in[2 ];
      sum += (coefs[6 ]*x + coefs[7 ]) * in[3 ];

      sum += (coefs[8 ]*x + coefs[9 ]) * in[4 ];
      sum += (coefs[10]*x + coefs[11]) * in[5 ];
      sum += (coefs[12]*x + coefs[13]) * in[6 ];
      sum += (coefs[14]*x + coefs[15]) * in[7 ];

      sum += (coefs[16]*x + coefs[17]) * in[8 ];
      sum += (coefs[18]*x + coefs[19]) * in[9 ];
      sum += (coefs[20]*x + coefs[21]) * in[10];
      sum += (coefs[22]*x + coefs[23]) * in[11];

      output[i] = sum;
    }
    fifo_read(&p->fifo, p->at.parts.integer, NULL);
    p->at.parts.integer = 0;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
}

static void u100_2(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t const * const input = (sox_sample_t *)fifo_read_ptr(&p->fifo) + p->pre;
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t * const output = fifo_reserve(output_fifo, max_num_out);
  sox_sample_t const * const coef0 = p->shared->poly_fir_coefs;
  {
    for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
      sox_sample_t const * const in = input + p->at.parts.integer;
      uint32_t fraction = p->at.parts.fraction;
      int phase = fraction >> (32 - u100_2_b); /* high-order bits */
      sox_sample_t const * const coefs = coef0 + 3*u100_l*phase;

      sox_sample_t x = (sox_sample_t) (fraction << u100_2_b) * (1 / MULT32);
      sox_sample_t sum = 0;

      sum += ((coefs[0 ]*x + coefs[1 ])*x + coefs[2 ]) * in[0 ];
      sum += ((coefs[3 ]*x + coefs[4 ])*x + coefs[5 ]) * in[1 ];
      sum += ((coefs[6 ]*x + coefs[7 ])*x + coefs[8 ]) * in[2 ];
      sum += ((coefs[9 ]*x + coefs[10])*x + coefs[11]) * in[3 ];

      sum += ((coefs[12]*x + coefs[13])*x + coefs[14]) * in[4 ];
      sum += ((coefs[15]*x + coefs[16])*x + coefs[17]) * in[5 ];
      sum += ((coefs[18]*x + coefs[19])*x + coefs[20]) * in[6 ];
      sum += ((coefs[21]*x + coefs[22])*x + coefs[23]) * in[7 ];

      sum += ((coefs[24]*x + coefs[25])*x + coefs[26]) * in[8 ];
      sum += ((coefs[27]*x + coefs[28])*x + coefs[29]) * in[9 ];
      sum += ((coefs[30]*x + coefs[31])*x + coefs[32]) * in[10];
      sum += ((coefs[33]*x + coefs[34])*x + coefs[35]) * in[11];

      output[i] = sum;
    }
    fifo_read(&p->fifo, p->at.parts.integer, NULL);
    p->at.parts.integer = 0;
  }
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
