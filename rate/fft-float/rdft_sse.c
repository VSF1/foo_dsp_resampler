/*
 * (I)RDFT transforms
 * Copyright (c) 2009 Alex Converse <alex dot converse at gmail dot com>
 * Copyright (c) 2018 lvqcl (intrinsics code)
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "ffti.h"
#include <xmmintrin.h>
#define PM128(x) (*(__m128*)(x))

static _MM_ALIGN16 const unsigned long N13[4] = { 0x00000000, 0x80000000, 0x00000000, 0x80000000 };
static _MM_ALIGN16 const float P1234[4] = { +0.5f, +0.5f, +0.5f, +0.5f };
static _MM_ALIGN16 const float P02M13[4] = { +0.5f, -0.5f, +0.5f, -0.5f };
static _MM_ALIGN16 const float P13M02[4] = { -0.5f, +0.5f, -0.5f, +0.5f };

#pragma warning(disable:4700)
void ff_rdft_calc_sse(RDFTContext *s, FFTSample *data, FFTComplex *tmp_buf)
{
    int i, i1, i2;
    const int n = 1 << s->nbits;
    const int n4 = n>>2;
    const float k1 = 0.5f;
    const float k2 = s->inverse ? -0.5f : 0.5f;
    const FFTSample *tcos = s->tcos;
    const FFTSample *tsin = s->tsin;

    const __m128 kk1 = PM128(P1234);
    const __m128 kk2 = s->inverse ? PM128(P02M13) : PM128(P13M02);
    const __m128 neg13 = PM128(N13);


    if (!s->inverse) {
        ff_fft_permute_sse(&s->fft, (FFTComplex*)data, tmp_buf);
        ff_fft_calc_sse(&s->fft, (FFTComplex*)data);
    }
    {
        FFTComplex ev, od, odsum;

        /* i=0 is a special case because of packing, the DC term is real, so we
           are going to throw the N/2 term (also real) in with it. */
        ev.re = data[0];
        data[0] = ev.re + data[1];
        data[1] = ev.re - data[1];

        ev.re =  k1 * (data[2] + data[n-2]);
        ev.im =  k1 * (data[3] - data[n-1]);
        od.im = -k2 * (data[2] - data[n-2]);
        od.re =  k2 * (data[3] + data[n-1]);
        if (s->negative_sin) {
            odsum.re = od.re * tcos[1] + od.im * tsin[1];
            odsum.im = od.im * tcos[1] - od.re * tsin[1];
        } else {
            odsum.re = od.re * tcos[1] - od.im * tsin[1];
            odsum.im = od.im * tcos[1] + od.re * tsin[1];
        }
        data[2  ] =  ev.re + odsum.re;
        data[3  ] =  ev.im + odsum.im;
        data[n-2] =  ev.re - odsum.re;
        data[n-1] = -ev.im + odsum.im;
    }

#if 0
    for (i=2, i1=4, i2=n-4; i < n4; i+=2, i1+=4, i2-=4)
    {
        FFTComplex ev1, od1, ev2, od2;

        ev1.re =  k1 * (data[i1  ] + data[i2  ]);
        ev1.im =  k1 * (data[i1+1] - data[i2+1]);
        ev2.re =  k1 * (data[i1+2] + data[i2-2]);
        ev2.im =  k1 * (data[i1+3] - data[i2-1]);

        od1.im = -k2 * (data[i1  ] - data[i2  ]);
        od1.re =  k2 * (data[i1+1] + data[i2+1]);
        od2.im = -k2 * (data[i1+2] - data[i2-2]);
        od2.re =  k2 * (data[i1+3] + data[i2-1]);

        if (s->negative_sin)
        {
            data[i2  ] = +(ev1.re - od1.re*tcos[i  ]) - od1.im*tsin[i  ];
            data[i2+1] = -(ev1.im - od1.im*tcos[i  ]) - od1.re*tsin[i  ];
            data[i2-2] = +(ev2.re - od2.re*tcos[i+1]) - od2.im*tsin[i+1];
            data[i2-1] = -(ev2.im - od2.im*tcos[i+1]) - od2.re*tsin[i+1];

            data[i1  ] =  ev1.re + od1.re*tcos[i  ] + od1.im*tsin[i  ];
            data[i1+1] =  ev1.im + od1.im*tcos[i  ] - od1.re*tsin[i  ];
            data[i1+2] =  ev2.re + od2.re*tcos[i+1] + od2.im*tsin[i+1];
            data[i1+3] =  ev2.im + od2.im*tcos[i+1] - od2.re*tsin[i+1];
        }
        else
        {
            data[i2  ] = +(ev1.re - od1.re*tcos[i  ]) + od1.im*tsin[i  ];
            data[i2+1] = -(ev1.im - od1.im*tcos[i  ]) + od1.re*tsin[i  ];
            data[i2-2] = +(ev2.re - od2.re*tcos[i+1]) + od2.im*tsin[i+1];
            data[i2-1] = -(ev2.im - od2.im*tcos[i+1]) + od2.re*tsin[i+1];

            data[i1  ] =  ev1.re + od1.re*tcos[i  ] - od1.im*tsin[i  ];
            data[i1+1] =  ev1.im + od1.im*tcos[i  ] + od1.re*tsin[i  ];
            data[i1+2] =  ev2.re + od2.re*tcos[i+1] - od2.im*tsin[i+1];
            data[i1+3] =  ev2.im + od2.im*tcos[i+1] + od2.re*tsin[i+1];
        }
    }
#else
    for (i=2, i1=4, i2=n-4; i < n4; i+=2, i1+=4, i2-=4)
    {
        __m128 dat0, dat1;
        __m128 ev, od, oq;
        __m128 mcos, msin;

        dat0 = _mm_load_ps(data+i1);                                 // d[i1]    d[i1+1]   d[i1+2]   d[i1+3]
        dat1 = _mm_loadu_ps(data+i2-2);                              // d[i2-2]  d[i2-1]   d[i2]     d[i2+1]
        dat1 = _mm_shuffle_ps(dat1, dat1, _MM_SHUFFLE(1,0,3,2));     // d[i2]    d[i2+1]   d[i2-2]   d[i2-1]
        dat1 = _mm_xor_ps(dat1, neg13);                              // d[i2]   -d[i2+1]   d[i2-2]  -d[i2-1]

        ev = _mm_add_ps(dat0, dat1);
        od = _mm_sub_ps(dat0, dat1);                        //    3        2        1        0
        ev = _mm_mul_ps(ev, kk1);                           // ev1.re   ev1.im   ev2.re   ev2.im
        od = _mm_mul_ps(od, kk2);                           // od1.im   od1.re   od2.im   od2.re
        oq = _mm_shuffle_ps(od, od, _MM_SHUFFLE(2,3,0,1));  // od1.re   od1.im   od2.re   od2.im

#ifdef _DEBUG
        mcos = msin = _mm_setzero_ps();
#endif
#pragma warning(disable:4700)
        mcos = _mm_loadl_pi(mcos, (__m64*)(tcos+i));
        msin = _mm_loadl_pi(msin, (__m64*)(tsin+i));
#pragma warning(default:4700)
        mcos = _mm_unpacklo_ps(mcos, mcos);
        msin = _mm_unpacklo_ps(msin, msin);
        oq = _mm_mul_ps(oq, mcos);
        od = _mm_mul_ps(od, msin);

        if (s->negative_sin)
        {
            dat1 = _mm_sub_ps(ev, oq);
            dat1 = _mm_xor_ps(dat1, neg13);
            dat1 = _mm_sub_ps(dat1, od);
            dat1 = _mm_shuffle_ps(dat1, dat1, _MM_SHUFFLE(1,0,3,2));

            od = _mm_xor_ps(od, neg13);
            dat0 = _mm_add_ps(ev, oq);
            dat0 = _mm_add_ps(dat0, od);

            _mm_storeu_ps(data+i2-2, dat1);
            _mm_store_ps(data+i1, dat0);
        }
        else
        {
            dat1 = _mm_sub_ps(ev, oq);
            dat1 = _mm_xor_ps(dat1, neg13);
            dat1 = _mm_add_ps(dat1, od);
            dat1 = _mm_shuffle_ps(dat1, dat1, _MM_SHUFFLE(1,0,3,2));

            od = _mm_xor_ps(od, neg13);
            dat0 = _mm_add_ps(ev, oq);
            dat0 = _mm_sub_ps(dat0, od);

            _mm_storeu_ps(data+i2-2, dat1);
            _mm_store_ps(data+i1, dat0);
        }
    }
#endif

    data[2*i+1] *= s->sign_convention;
    if (s->inverse) {
        data[0] *= k1;
        data[1] *= k1;
        ff_fft_permute_sse(&s->fft, (FFTComplex*)data, tmp_buf);
        ff_fft_calc_sse(&s->fft, (FFTComplex*)data);
    }
}
