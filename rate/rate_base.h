/* Effect: change sample rate  Copyright (c) 2008,12,16 robs@users.sourceforge.net
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

/* Inspired by, and builds upon some of the ideas presented in:
 * `The Quest For The Perfect Resampler' by Laurent De Soras;
 * http://ldesoras.free.fr/doc/articles/resampler-en.pdf
 */

#ifdef RATE_SINGLE_PRECISION

#ifdef RATE_DOUBLE_PRECISION
#error select single or double precision
#endif

#ifdef SSE3_
#error single precision uses only SSE1 instructions so it works (almost) everywhere
#endif

#endif

#include <assert.h>
#include <string.h>

#include "sox_i.h"

#define  FIFO_SIZE_T int
#include "fifo.h"

#ifdef RATE_SINGLE_PRECISION

#ifdef SSE_
#define lsx_safe_rdft ff_rdft_SSE
#else
#define lsx_safe_rdft ff_rdft_generic
#endif

#else

#ifdef SSE3_
#define lsx_safe_rdft lsx_safe_rdft_SSE3
#else
#define lsx_safe_rdft lsx_safe_rdft_generic
#endif

#endif

#define raw_coef_t double

#ifdef RATE_SINGLE_PRECISION
  #define sox_sample_t float
#else
  #define sox_sample_t double
#endif

/* there may me a phase shift if num_coefs != num_coefs4; didn't look further yet */
#ifdef CREATE_4X_NUMTAPS
  #define num_coefs4 num_coefs
#else
  #define num_coefs4 ((num_coefs + 3) & ~3) /* align coefs for SSE */
#endif

#if defined M_PIl
  #define hi_prec_clock_t long double /* __float128 is also a (slow) option */
#else
  #define hi_prec_clock_t double
#endif


typedef struct {
  int dft_length, num_taps, post_peak;
  sox_sample_t * coefs;
  FFTComplex * tmp_buf; /* for ffmpeg FFT permute */
} dft_filter_t;

typedef struct { /* So generated filter coefs may be shared between channels */
  sox_sample_t  * poly_fir_coefs;
  dft_filter_t dft_filter[2];
} rate_shared_t;

struct stage;
typedef void (* stage_fn_t)(struct stage * input, fifo_t * output);
typedef struct stage {
  /* Common to all stage types: */
  stage_fn_t fn;
  fifo_t     fifo;
  int        pre;       /* Number of past samples to store */
  int        pre_post;  /* pre + number of future samples to store */
  int        preload;   /* Number of zero samples to pre-load the fifo */
  double     out_in_ratio; /* For buffer management. */

  /* For a stage with variable (run-time generated) filter coefs: */
  rate_shared_t * shared;
  int        dft_filter_num; /* Which, if any, of the 2 DFT filters to use */

  /* For a stage with variable L/M: */
  union {               /* 32bit.32bit fixed point arithmetic */
    #if defined(WORDS_BIGENDIAN)
    struct {int32_t integer; uint32_t fraction;} parts;
    #else
    struct {uint32_t fraction; int32_t integer;} parts;
    #endif
    int64_t all;
#ifdef RATE_SINGLE_PRECISION
    #define MULT32 (65536.f * 65536.f)
#else
    #define MULT32 (65536. * 65536.)
#endif

    hi_prec_clock_t hi_prec_clock;
  } at, step;
  sox_bool   use_hi_prec_clock;
  int        L, remL, remM;
  int        n, phase_bits;
} stage_t;

#define stage_occupancy(s) max(0, fifo_occupancy(&(s)->fifo) - (s)->pre_post)
#define stage_read_p(s) ((sox_sample_t *)fifo_read_ptr(&(s)->fifo) + (s)->pre)

#if 0
static void cubic_stage_fn(stage_t * p, fifo_t * output_fifo)
{
  int i, num_in = stage_occupancy(p), max_num_out = 1 + num_in*p->out_in_ratio;
  sox_sample_t const * input = stage_read_p(p);
  sox_sample_t * output = fifo_reserve(output_fifo, max_num_out);

  for (i = 0; p->at.parts.integer < num_in; ++i, p->at.all += p->step.all) {
    sox_sample_t const * s = input + p->at.parts.integer;
    sox_sample_t x = p->at.parts.fraction * (1 / MULT32);
    sox_sample_t b = .5*(s[1]+s[-1])-*s, a = (1/6.)*(s[2]-s[1]+s[-1]-*s-4*b);
    sox_sample_t c = s[1]-*s-a-b;
    output[i] = ((a*x + b)*x + c)*x + *s;
  }
  assert(max_num_out - i >= 0);
  fifo_trim_by(output_fifo, max_num_out - i);
  fifo_read(&p->fifo, p->at.parts.integer, NULL);
  p->at.parts.integer = 0;
}
#endif

static void dft_stage_fn(stage_t * p, fifo_t * output_fifo);

static void dft_stage_init(
    unsigned instance, double Fp, double Fs, double Fn, double att,
    double phase, stage_t * stage, int L, int M)
{
  dft_filter_t * f = &stage->shared->dft_filter[instance];

  if (!f->num_taps) {
    int num_taps = 0, dft_length, i;
    int k = phase == 50 && lsx_is_power_of_2(L) && Fn == L? L << 1 : 4;
    double * h = lsx_design_lpf(Fp, Fs, Fn, att, &num_taps, -k, -1.);

    if (phase != 50)
      lsx_fir_to_phase(&h, &num_taps, &f->post_peak, phase);
    else f->post_peak = num_taps / 2;

    dft_length = lsx_set_dft_length(num_taps);
    f->coefs = lsx_aligned_calloc(dft_length, sizeof(*f->coefs));
    for (i = 0; i < num_taps; ++i)
      f->coefs[(i + dft_length - num_taps + 1) & (dft_length - 1)]
        = h[i] / dft_length * 2 * L;
    lsx_free(h);
    f->num_taps = num_taps;
    f->dft_length = dft_length;
#ifdef RATE_SINGLE_PRECISION
    f->tmp_buf = lsx_aligned_malloc(dft_length*sizeof(FFTComplex)/2); // for ffmpeg FFT out-of-place permute routine
#else
    f->tmp_buf = NULL;
#endif
    lsx_safe_rdft(dft_length, 1, f->coefs, f->tmp_buf);
  }
  stage->fn = dft_stage_fn;
  stage->preload = f->post_peak / L;
  stage->remL    = f->post_peak % L;
  stage->L = L;
  stage->step.parts.integer = abs(3-M) == 1 && Fs == 1? -M/2 : M;
  stage->dft_filter_num = instance;
}

static __inline div_t __cdecl div1(int numer, int denom)
{
    div_t result;

    result.quot = numer / denom;
    result.rem = numer % denom;

    return result;
}

#ifdef SSE_
#include <xmmintrin.h>
#endif
#ifdef SSE3_
#include <pmmintrin.h>
#endif
#include <stdint.h>

#include "dft_filter.h"
#include "prepare_coefs.h"

#if defined(SSE3_)
#include "rate_filters_SSE3.h"
#elif defined(SSE_)
#include "rate_filters_SSE.h"
#else
#include "rate_filters_generic.h"
#endif


typedef struct {
  double     factor;
  size_t     samples_in, samples_out;
  int        num_stages;
  stage_t    * stages;

  size_t in_samplerate, out_samplerate;
} rate_t;

#define pre_stage       p->stages[shift]
#define arb_stage       p->stages[shift + have_pre_stage]
#define post_stage      p->stages[shift + have_pre_stage + have_arb_stage]
#define have_pre_stage  (preM  * preL  != 1)
#define have_arb_stage  (arbM  * arbL  != 1)
#define have_post_stage (postM * postL != 1)

#define TO_3dB(a)       ((1.6e-6*a-7.5e-4)*a+.646)
#define LOW_Q_BW0_PC    (67 + 5 / 8.)

typedef enum {
  rolloff_none, rolloff_small /* <= 0.01 dB */, rolloff_medium /* <= 0.35 dB */
} rolloff_t;

static void rate_init(
  /* Private work areas (to be supplied by the client):                       */
  rate_t * p,                /* Per audio channel.                            */
  rate_shared_t * shared,    /* Between channels (undergoing same rate change)*/

  /* Public parameters:                                             Typically */
  double factor,             /* Input rate divided by output rate.            */
  double bits,               /* Required bit-accuracy (pass + stop)  16|20|28 */
  double phase,              /* Linear/minimum etc. filter phase.       50    */
  double bw_pc,              /* Pass-band % (0dB pt.) to preserve.   91.3|98.4*/
  double anti_aliasing_pc,   /* % bandwidth without aliasing            100   */
  rolloff_t rolloff,         /* Pass-band roll-off                    small   */
  sox_bool maintain_3dB_pt,  /*                                        true   */

  /* Primarily for test/development purposes:                                 */
  sox_bool use_hi_prec_clock,/* Increase irrational ratio accuracy.   false   */
  int interpolator,          /* Force a particular coef interpolator.   -1    */
  int max_coefs_size,        /* k bytes of coefs to try to keep below.  400   */
  sox_bool noSmallIntOpt)    /* Disable small integer optimisations.  false   */
{
  double att = (bits + 1) * linear_to_dB(2.), attArb = att;    /* pass + stop */
  double tbw0 = 1 - bw_pc / 100, Fs_a = 2 - anti_aliasing_pc / 100;
  double arbM = factor, tbw_tighten = 1;
  int n = 0, i, preL = 1, preM = 1, shift = 0, arbL = 1, postL = 1, postM = 1;
  sox_bool upsample = sox_false, rational = sox_false, iOpt = !noSmallIntOpt;
  int mode = rolloff > rolloff_small? factor > 1 || bw_pc > LOW_Q_BW0_PC :
    ceil(2 + (bits - 17) / 4);
  stage_t * s;

  assert(factor > 0);
  assert(/*!bits || */(15 <= bits && bits <= 33));
  assert(0 <= phase && phase <= 100);
  assert(53 <= bw_pc && bw_pc <= 100);
  assert(85 <= anti_aliasing_pc && anti_aliasing_pc <= 100);

  p->factor = factor;
  /*if (bits)*/ while (!n++) {                               /* Determine stages: */
    int try_i, L, M, x, maxL = interpolator > 0? 1 : mode? 2048 :
      ceil(max_coefs_size * 1000. / (U100_l * sizeof(sox_sample_t)));
    double d, epsilon = 0, frac;
    upsample = arbM < 1;
    for (i = arbM * .5, shift = 0; i >>= 1; arbM *= .5, ++shift);
    preM = upsample || (arbM > 1.5 && arbM < 2);
    postM = 1 + (arbM > 1 && preM), arbM /= postM;
    preL = 1 + (!preM && arbM < 2) + (upsample && mode), arbM *= preL;
    if ((frac = arbM - (int)arbM))
      epsilon = fabs(floor(frac * MULT32 + .5) / (frac * MULT32) - 1);
    for (i = 1, rational = !frac; i <= maxL && !rational; ++i) {
      d = frac * i, try_i = d + .5;
      if ((rational = fabs(try_i / d - 1) <= epsilon)) {    /* No long doubles! */
        if (try_i == i)
          arbM = ceil(arbM), shift += x = arbM > 3, arbM /= 1 + x;
        else arbM = i * (int)arbM + try_i, arbL = i;
      }
    }
    L = preL * arbL, M = arbM * postM, x = (L|M)&1, L >>= !x, M >>= !x;
    if (iOpt && postL == 1 && (d = preL * arbL / arbM) > 4 && d != 5) {
      for (postL = 4, i = d / 16; i >>= 1; postL <<= 1);
      arbM = arbM * postL / arbL / preL, arbL = 1, n = 0;
    } else if (rational && (max(L, M) < 3 + 2 * iOpt || L * M < 6 * iOpt))
      preL = L, preM = M, arbM = arbL = postM = 1;
    if (!mode && (!rational || !n))
      ++mode, n = 0;
  }

  p->num_stages = shift + have_pre_stage + have_arb_stage + have_post_stage;
  p->stages = lsx_calloc(p->num_stages + 1, sizeof(*p->stages));
  for (i = 0; i < p->num_stages; ++i)
    p->stages[i].shared = shared;

  if ((n = p->num_stages) > 1) {                              /* Att. budget: */
    if (have_arb_stage)
      att += linear_to_dB(2.), attArb = att, --n;
    att += linear_to_dB((double)n);
  }

  for (n = 0; n + 1u < array_length(half_firs) && att > half_firs[n].att; ++n);
  for (i = 0, s = p->stages; i < shift; ++i, ++s) {
    s->fn = half_firs[n].fn;
    s->pre_post = 4 * half_firs[n].num_coefs;
    s->preload = s->pre = s->pre_post >> 1;
  }

  if (have_pre_stage) {
    if (maintain_3dB_pt && have_post_stage) {    /* Trans. bands overlapping. */
      double tbw3 = tbw0 * TO_3dB(att);               /* TODO: consider Fs_a. */
      double x = ((2.1429e-4 - 5.2083e-7 * att) * att - .015863) * att + 3.95;
      x = att * pow((tbw0 - tbw3) / (postM / (factor * postL) - 1 + tbw0), x);
      if (x > .035) {
        tbw_tighten = ((4.3074e-3 - 3.9121e-4 * x) * x - .040009) * x + 1.0014;
      }
    }
    dft_stage_init(0, 1 - tbw0 * tbw_tighten, Fs_a, preM? max(preL, preM) :
        arbM / arbL, att, phase, &pre_stage, preL, max(preM, 1));
  }

  /*if (!bits) {                                   Quick and dirty arb stage:
    arb_stage.fn = cubic_stage_fn;
    arb_stage.step.all = arbM * MULT32 + .5;
    arb_stage.pre_post = max(3, arb_stage.step.parts.integer);
    arb_stage.preload = arb_stage.pre = 1;
    arb_stage.out_in_ratio = MULT32 * arbL / arb_stage.step.all;
  }
  else*/ if (have_arb_stage) {                     /* Higher quality arb stage: */
    poly_fir_t const * f = &poly_firs[6*(upsample + !!preM) + mode - !upsample];
    int order, num_coefs = f->interp[0].scalar, phase_bits, phases, coefs_size;
    double x = .5, at, Fp, Fs, Fn, mult = upsample? 1 : arbL / arbM;
    poly_fir1_t const * f1;

    Fn = !upsample && preM? x = arbM / arbL : 1;
    Fp = !preM? mult : mode? .5 : 1;
    Fs = 2 - Fp;           /* Ignore Fs_a; it would have little benefit here. */
    Fp *= 1 - tbw0;
    if (rolloff > rolloff_small && mode)
      Fp = !preM? mult * .5 - .125 : mult * .05 + .1;
    else if (rolloff == rolloff_small)
      Fp = Fs - (Fs - .148 * x - Fp * .852) * (.00813 * bits + .973);

    i = (interpolator < 0? !rational : max(interpolator, !rational)) - 1;
    do {
      f1 = &f->interp[++i];
      assert(f1->fn);
      if (i)
        arbM /= arbL, arbL = 1, rational = sox_false;
      phase_bits = ceil(f1->scalar + log(mult)/log(2.));
      phases = !rational? (1 << phase_bits) : arbL;
      if (!f->interp[0].scalar) {
        int phases0 = max(phases, 19), n0 = 0;
        lsx_design_lpf(Fp, Fs, -Fn, attArb, &n0, phases0, f->beta);
        num_coefs = n0 / phases0 + 1, num_coefs += num_coefs & !preM;
      }
      if ((num_coefs & 1) && rational && (arbL & 1))
        phases <<= 1, arbL <<= 1, arbM *= 2;
      at = arbL * .5 * (num_coefs & 1);
      order = i + (i && mode > 4);
      coefs_size = num_coefs4 * phases * (order + 1) * sizeof(sox_sample_t);
    } while (interpolator < 0 && i < 2 && f->interp[i+1].fn &&
        coefs_size / 1000 > max_coefs_size);

    if (!arb_stage.shared->poly_fir_coefs) {
      int num_taps = num_coefs * phases - 1;
      raw_coef_t * coefs = lsx_design_lpf(
          Fp, Fs, Fn, attArb, &num_taps, phases, f->beta);
      arb_stage.shared->poly_fir_coefs = prepare_coefs(
          coefs, num_coefs, phases, order, 1);
      lsx_free(coefs);
    }
    arb_stage.fn = f1->fn;
    arb_stage.pre_post = num_coefs4 - 1;
    arb_stage.preload = (num_coefs - 1) >> 1;
    arb_stage.n = num_coefs4; assert(arb_stage.n%4 == 0);
    arb_stage.phase_bits = phase_bits;
    arb_stage.L = arbL;
    arb_stage.use_hi_prec_clock = mode > 1 && use_hi_prec_clock && !rational;
    if (arb_stage.use_hi_prec_clock) {
      arb_stage.at.hi_prec_clock = at;
      arb_stage.step.hi_prec_clock = arbM;
      arb_stage.out_in_ratio = arbL / arb_stage.step.hi_prec_clock;
    } else {
      arb_stage.at.all = at * MULT32 + .5;
      arb_stage.step.all = arbM * MULT32 + .5;
      arb_stage.out_in_ratio = MULT32 * arbL / arb_stage.step.all;
    }
  }

  if (have_post_stage)
    dft_stage_init(1, 1 - (1 - (1 - tbw0) *
        (upsample? factor * postL / postM : 1)) * tbw_tighten, Fs_a,
        (double)max(postL, postM), att, phase, &post_stage, postL, postM);

  for (i = 0, s = p->stages; i < p->num_stages; ++i, ++s) {
    fifo_create(&s->fifo, (int)sizeof(sox_sample_t));
    //memset(fifo_reserve(&s->fifo, s->preload), 0, sizeof(sox_sample_t)*s->preload);
    fifo_write0(&s->fifo, s->preload);
  }
  fifo_create(&s->fifo, (int)sizeof(sox_sample_t));
}

static void rate_process(rate_t * p)
{
  stage_t * stage = p->stages;
  int i;

  for (i = 0; i < p->num_stages; ++i, ++stage)
    stage->fn(stage, &(stage+1)->fifo);
}

static sox_sample_t * rate_input(rate_t * p, sox_sample_t const * samples, size_t n)
{
  p->samples_in += n;
  while ((p->samples_in > p->in_samplerate) && (p->samples_out > p->out_samplerate))
  {
      p->samples_in -= p->in_samplerate;
      p->samples_out -= p->out_samplerate;
  }
  return fifo_write(&p->stages[0].fifo, (int)n, samples);
}

static sox_sample_t const * rate_output(rate_t * p, sox_sample_t * samples, size_t * n)
{
  fifo_t * fifo = &p->stages[p->num_stages].fifo;
  p->samples_out += *n = min(*n, (size_t)fifo_occupancy(fifo));
  return fifo_read(fifo, (int)*n, samples);
}

static sox_sample_t buff[1024] = {0};

static void rate_flush(rate_t * p)
{
  fifo_t * fifo = &p->stages[p->num_stages].fifo;
  size_t samples_out = (size_t)(p->samples_in / p->factor + .5);

  if (samples_out > p->samples_out) {
    size_t remaining = samples_out - p->samples_out;
    while ((size_t)fifo_occupancy(fifo) < remaining) {
      rate_input(p, buff, (size_t) 1024);
      rate_process(p);
    }
    fifo_trim_to(fifo, (int)remaining);
    p->samples_out = p->samples_in = 0;
  }
}

static void rate_close_shared(rate_shared_t * shared)
{
  lsx_aligned_free(shared->dft_filter[0].coefs);
  lsx_aligned_free(shared->dft_filter[1].coefs);
  lsx_aligned_free(shared->poly_fir_coefs);

  lsx_aligned_free(shared->dft_filter[0].tmp_buf);
  lsx_aligned_free(shared->dft_filter[1].tmp_buf);

  memset(shared, 0, sizeof(*shared));
}

static void rate_close(rate_t * p)
{
  int i;

  for (i = 0; i <= p->num_stages; ++i)
    fifo_delete(&p->stages[i].fifo);

  lsx_free(p->stages);
}

/*------------------------------- Wrapper for foobar2000 --------------------------------*
 *                                                         Copyright (c) 2008-2012 lvqcl */

#include "rate_i.h"

typedef struct RR_internal_tag
{
    RR_vtable x;

    int nchannels;
    rate_t* rate;
    rate_shared_t shared;
    size_t isamp_max;
} RR_internal;


typedef struct priv_tag {
  int  coef_interp, max_coefs_size;
  rolloff_t rolloff;
  double  bit_depth, phase, bw_0dB_pc, anti_aliasing_pc;
  sox_bool  use_hi_prec_clock, noIOpt, given_0dB_pt;
} priv;

static void convert_settings(const RR_config* cfg, priv* w);

static int RR_init_x(RR_handle* h, const RR_config* config, int nchannels)
{
    RR_internal* p = (RR_internal*)h;
    priv w;
    int i;
    double factor;
    if (p == NULL) return RR_NULLHANDLE;

    convert_settings(config, &w);

    factor = (double)config->in_rate/(double)config->out_rate;
    if (factor > 5644.8 || factor < 1.0/5644.8) return RR_INVPARAM;
    p->rate = (rate_t*) lsx_calloc(nchannels, sizeof(rate_t));
    p->nchannels = nchannels;
    p->isamp_max = 1048576; if (factor < 1) p->isamp_max *= factor; /*if (p->isamp_max < 1) p->isamp_max = 1;*/

    for(i=0; i < nchannels; i++)
    {
        p->rate[i].in_samplerate = config->in_rate;
        p->rate[i].out_samplerate = config->out_rate;
        rate_init(&p->rate[i], &p->shared, factor,
            w.bit_depth, w.phase, w.bw_0dB_pc, w.anti_aliasing_pc, w.rolloff, !w.given_0dB_pt,
            w.use_hi_prec_clock, w.coef_interp, w.max_coefs_size, w.noIOpt);
    }

    return RR_OK;
}

static void RR_close_x(RR_handle* h)
{
    RR_internal* p = (RR_internal*)h;
    int i;
    if (p == NULL) return;

    for(i=0; i < p->nchannels; i++)
        rate_close(&p->rate[i]);
    lsx_free(p->rate);
    rate_close_shared(&p->shared);

    lsx_free(p);
}

static __inline void interleave(fb_sample_t* out, const sox_sample_t* in, size_t nch, size_t nsamples)
{
    size_t n;
    for (n = 0; n < nsamples; ++n) out[n*nch] = (fb_sample_t) in[n];
}

static __inline void deinterleave(sox_sample_t* out, const fb_sample_t* in, size_t nch, size_t nsamples)
{
    size_t n;
    for (n = 0; n < nsamples; ++n) out[n] = (sox_sample_t) in[n*nch];
}

static int RR_flow_x(RR_handle *h, const fb_sample_t* const ibuf, fb_sample_t* const obuf, size_t isamp, size_t osamp, size_t* i_used, size_t* o_done)
{
    RR_internal* p = (RR_internal*)h;
    size_t itmpdummy;
    int i;
    int nchan;
    sox_sample_t const * s;
    size_t odone, odone2;
    assert(p);
    if (p == NULL) return RR_NULLHANDLE;

    if (i_used == NULL) i_used = &itmpdummy;
    if (ibuf == NULL) { isamp = 0; *i_used = 0; }
    if (isamp > p->isamp_max) isamp = p->isamp_max;
    nchan = p->nchannels;

    for(i=0; i < nchan; i++)
    {
        fb_sample_t* const obuf_i = obuf + i;
        *i_used = 0;
        odone = osamp;
        s = rate_output(&p->rate[i], NULL, &odone);
        interleave(obuf_i, s, nchan, odone);
        *o_done = odone;
        if (isamp)
        {
            const fb_sample_t* const ibuf_i = ibuf + i;
            sox_sample_t * t = rate_input(&p->rate[i], NULL, isamp);
            deinterleave(t, ibuf_i, nchan, isamp);
            rate_process(&p->rate[i]);
            *i_used = isamp;
        }
        if (odone < osamp)
        {
            fb_sample_t* const obuf_2 = obuf_i + odone*nchan;
            odone2 = osamp - odone;
            s = rate_output(&p->rate[i], NULL, &odone2);
            interleave(obuf_2, s, nchan, odone2);
            *o_done += odone2;
        }
    }

    return RR_OK;
}

static int RR_push_x(RR_handle *h, const fb_sample_t* const ibuf, size_t isamp)
{
    RR_internal* p = (RR_internal*)h;
    int i, nchan;
    assert(p);
    if (p == NULL) return RR_NULLHANDLE;

    if (ibuf == NULL || isamp == 0) return RR_OK;
    if (isamp > p->isamp_max) isamp = p->isamp_max;
    nchan = p->nchannels;

    for(i=0; i < nchan; i++)
    {
        const fb_sample_t* const ibuf_i = ibuf + i;
        sox_sample_t * t = rate_input(&p->rate[i], NULL, isamp);
        deinterleave(t, ibuf_i, nchan, isamp);
        rate_process(&p->rate[i]);
    }

    return RR_OK;
}

static int RR_pull_x(RR_handle *h, fb_sample_t* const obuf, size_t osamp, size_t* o_done)
{
    RR_internal* p = (RR_internal*)h;
    int i, nchan;
    sox_sample_t const * s;
    size_t odone;
    assert(p);
    if (p == NULL) return RR_NULLHANDLE;

    if (obuf == NULL || osamp == 0) { if (o_done) *o_done = 0; return RR_OK; }
    nchan = p->nchannels;

    for(i=0; i < nchan; i++)
    {
        fb_sample_t* const obuf_i = obuf + i;
        odone = osamp;
        s = rate_output(&p->rate[i], NULL, &odone);
        interleave(obuf_i, s, nchan, odone);
    }
    if (o_done) *o_done = odone;

    return RR_OK;
}

static int RR_drain_x(RR_handle *h)
{
    RR_internal* p = (RR_internal*)h;
    int i;
    assert(p);
    if (p == NULL) return RR_NULLHANDLE;

    for(i=0; i < p->nchannels; i++)
        rate_flush(&p->rate[i]);
    return RR_OK;
}

static void convert_settings(const RR_config* cfg, priv* p)
{
    double rej, bw_3dB_pc;
    int quality;
    sox_bool allow_aliasing;

    memset(p, 0, sizeof(*p));

    bw_3dB_pc = cfg->bandwidth;
    allow_aliasing = (sox_bool)cfg->allow_aliasing;
    p->phase = cfg->phase;
    p->coef_interp = -1;
    p->max_coefs_size = 400;

    if (cfg->quality == RR_best) {
        quality = 6;
        p->rolloff = rolloff_none;
    }
    else {
        quality = 4;
        p->rolloff = rolloff_small;
    }
    p->anti_aliasing_pc = allow_aliasing? bw_3dB_pc : 100;

    p->bit_depth = 16 + 4 * max(quality - 3, 0);
    rej = p->bit_depth * linear_to_dB(2.);
    p->bw_0dB_pc = 100 - (100 - bw_3dB_pc) / TO_3dB(rej);
    //p->bw_0dB_pc = quality == 1 ? LOW_Q_BW0_PC : p->bw_0dB_pc;

    p->use_hi_prec_clock = p->noIOpt = p->given_0dB_pt = sox_false;
}


#ifdef RATE_SINGLE_PRECISION

#ifdef SSE_
#define RR_ctor_ZZZ RR_ctor_SSE
#else
#define RR_ctor_ZZZ RR_ctor_float
#endif

#else

#ifdef SSE3_
#define RR_ctor_ZZZ RR_ctor_SSE3
#else
#define RR_ctor_ZZZ RR_ctor_double
#endif

#endif

RR_handle* RR_ctor_ZZZ(const RR_config* config, int nchannels)
{
    RR_internal* p;
    p = (RR_internal*)lsx_calloc(1, sizeof (RR_internal) );  //here *p.shared.dft_filter[n].coefs and .tmp_buf are initialized by NULL
    if (p == NULL) return NULL;

    p->x.init = RR_init_x;
    p->x.flow = RR_flow_x;
    p->x.push = RR_push_x;
    p->x.pull = RR_pull_x;
    p->x.drain = RR_drain_x;
    p->x.close = RR_close_x;

    RR_init_x((RR_handle*)p, config, nchannels);

    return (RR_handle*)p;
}
