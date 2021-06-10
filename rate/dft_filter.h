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

#ifdef SSE_

static _MM_ALIGN16 const uint32_t C_PNPN[4] = { 0x00000000, 0x80000000, 0x00000000, 0x80000000 };
#define SIGN _mm_load_ps((float*)C_PNPN)

static __inline __m128 ZMUL2(__m128 ai, __m128 b, __m128 sign)
{
    __m128 ar = ai;
    ar = _mm_shuffle_ps(ar, ar, _MM_SHUFFLE(2, 2, 0, 0));
    ai = _mm_shuffle_ps(ai, ai, _MM_SHUFFLE(3, 3, 1, 1));
    ar = _mm_mul_ps(ar, b);
    ai = _mm_mul_ps(ai, b);
    ai = _mm_xor_ps(ai, sign);
    ai = _mm_shuffle_ps(ai, ai, _MM_SHUFFLE(2, 3, 0, 1));
    return _mm_add_ps(ar, ai);
}

#endif

#ifdef SSE3_

#define _mm_shuffle_lh(x, y) (_mm_shuffle_pd((x), (y), _MM_SHUFFLE2(0, 1)))

static __inline __m128d ZMUL(__m128d ai, __m128d b)
{
    __m128d ar=ai, b2=b;

    ar = _mm_movedup_pd(ai);
    ai = _mm_unpackhi_pd(ai, ai);
    b2 = _mm_shuffle_lh(b2, b2);
    ar = _mm_mul_pd(ar, b);
    ai = _mm_mul_pd(ai, b2);

    return _mm_addsub_pd(ar, ai);
}

#endif



static void dft_stage_fn(stage_t * p, fifo_t * output_fifo)
{
  sox_sample_t * output;
#ifdef SSE_
  __m128 coef, outp;
  __m128 sign = SIGN;
#endif
#ifdef SSE3_
  __m128d coef, outp;
#else
  sox_sample_t tmp;
#endif
  int i, j, num_in = max(0, fifo_occupancy(&p->fifo));
  rate_shared_t const * s = p->shared;
  dft_filter_t const * f = &s->dft_filter[p->dft_filter_num];
  int const overlap = f->num_taps - 1;
  const sox_sample_t * const coeff = f->coefs;

  while (p->remL + p->L * num_in >= f->dft_length)
  {
    div_t divd = div1(f->dft_length - overlap - p->remL + p->L - 1, p->L);
    sox_sample_t const * input = fifo_read_ptr(&p->fifo);
    fifo_read(&p->fifo, divd.quot, NULL);
    num_in -= divd.quot;

    output = fifo_reserve_aligned(output_fifo, f->dft_length);
    if (lsx_is_power_of_2(p->L)) /* F-domain */
    {
      int portion = f->dft_length / p->L;
      memcpy(output, input, (unsigned)portion * sizeof(*output));
      lsx_safe_rdft(portion, 1, output, f->tmp_buf);
      for (i = portion + 2; i < (portion << 1); i += 2)
      {
        output[i  ] =  output[(portion << 1) - i];
        output[i+1] = -output[(portion << 1) - i + 1];
      }
      output[portion  ] = output[1];
      output[portion+1] = 0;
      output[1] = output[0];
      for (portion <<= 1; i < f->dft_length; i += portion, portion <<= 1)
      {
        memcpy(output + i, output, portion * sizeof(*output));
        output[i + 1] = 0;
      }
    }
    else
    {
      if (p->L == 1)
        memcpy(output, input, f->dft_length * sizeof(*output));
      else
      {
        memset(output, 0, f->dft_length * sizeof(*output));
        for (j = 0, i = p->remL; i < f->dft_length; ++j, i += p->L)
          output[i] = input[j];
        p->remL = p->L - 1 - divd.rem;
      }
      lsx_safe_rdft(f->dft_length, 1, output, f->tmp_buf);
    }
    output[0] *= coeff[0];
    if (p->step.parts.integer > 0)
    {
      output[1] *= coeff[1];
#if defined(SSE3_)
      for (i = 2; i < f->dft_length; i += 2)
      {
        outp = _mm_load_pd(output+i);
        coef = _mm_load_pd(coeff+i);
        _mm_store_pd(output+i, ZMUL(outp, coef));
      }
#elif defined(SSE_)
      tmp = output[2];
      output[2] = coeff[2] * tmp - coeff[3] * output[3];
      output[3] = coeff[3] * tmp + coeff[2] * output[3];
      for (i = 4; i < f->dft_length; i += 4)
      {
        outp = _mm_load_ps(output+i);
        coef = _mm_load_ps(coeff+i);
        _mm_store_ps(output+i, ZMUL2(outp, coef, sign));
      }
#else
      for (i = 2; i < f->dft_length; i += 2)
      {
        tmp = output[i];
        output[i  ] = coeff[i  ] * tmp - coeff[i+1] * output[i+1];
        output[i+1] = coeff[i+1] * tmp + coeff[i  ] * output[i+1];
      }
#endif
      lsx_safe_rdft(f->dft_length, -1, output, f->tmp_buf);
      if (p->step.parts.integer != 1)
      {
        for (j = 0, i = p->remM; i < f->dft_length - overlap; ++j, i += p->step.parts.integer)
          output[j] = output[i];
        p->remM = i - (f->dft_length - overlap);
        fifo_trim_by(output_fifo, f->dft_length - j);
      }
      else fifo_trim_by(output_fifo, overlap);
    }
    else /* F-domain */
    {
      int m = -p->step.parts.integer;
#if defined(SSE3_)
      for (i = 2; i < (f->dft_length >> m); i += 2)
      {
        outp = _mm_load_pd(output+i);
        coef = _mm_load_pd(coeff+i);
        _mm_store_pd(output+i, ZMUL(outp, coef));
      }
#elif defined(SSE_)
      tmp = output[2];
      output[2] = coeff[2] * tmp - coeff[3] * output[3];
      output[3] = coeff[3] * tmp + coeff[2] * output[3];
      for (i = 4; i < (f->dft_length >> m); i += 4)
      {
        outp = _mm_load_ps(output+i);
        coef = _mm_load_ps(coeff+i);
        _mm_store_ps(output+i, ZMUL2(outp, coef, sign));
      }
#else
      for (i = 2; i < (f->dft_length >> m); i += 2)
      {
        tmp = output[i];
        output[i  ] = coeff[i  ] * tmp - coeff[i+1] * output[i+1];
        output[i+1] = coeff[i+1] * tmp + coeff[i  ] * output[i+1];
      }
#endif
      output[1] = coeff[i] * output[i] - coeff[i+1] * output[i+1];
      lsx_safe_rdft(f->dft_length >> m, -1, output, f->tmp_buf);
      fifo_trim_by(output_fifo, (((1 << m) - 1) * f->dft_length + overlap) >>m);
    }
  }
}
