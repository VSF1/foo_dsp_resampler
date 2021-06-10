/* Copyright(C) Takuya OOURA
 * email: ooura@mmm.t.u-tokyo.ac.jp
 * download: http://momonga.t.u-tokyo.ac.jp/~ooura/fft.html
 * You may use, copy, modify this code for any purpose and
 * without fee. You may distribute this ORIGINAL package.
 *
 *
 * + some changes from SoX project
 * + SSE accelerated functions: Copyright (C) lvqcl
*/

#include "fft4g.h"

#include <pmmintrin.h>
#include <stdint.h>

#define _mm_shuffle_lh(x, y) (_mm_shuffle_pd((x), (y), _MM_SHUFFLE2(0, 1)))

static _MM_ALIGN16 const uint64_t CC_NP[2] = { 0x8000000000000000ULL, 0x0000000000000000ULL };
static _MM_ALIGN16 const uint64_t CC_PN[2] = { 0x0000000000000000ULL, 0x8000000000000000ULL };

#define C_NP ((double*)CC_NP)
#define C_PN ((double*)CC_PN)

#define _mm_neg_lo(x) (_mm_xor_pd((x), *(__m128d*)(CC_NP)))
#define _mm_neg_hi(x) (_mm_xor_pd((x), *(__m128d*)(CC_PN)))

#define _mm_neg_lo2(x) (_mm_xor_pd((x), (neg_lo)))
#define _mm_neg_hi2(x) (_mm_xor_pd((x), (neg_hi)))


static __inline __m128d zmul(__m128d ai, __m128d b)
{
    __m128d ar = ai;

    ar = _mm_movedup_pd(ar);      /* ar = [a.r a.r] */
    ai = _mm_unpackhi_pd(ai, ai); /* ai = [a.i a.i] */
    ar = _mm_mul_pd(ar, b);       /* ar = [a.r*b.r a.r*b.i] */
    b  = _mm_shuffle_lh(b, b);    /* b = [b.i b.r] */
    ai = _mm_mul_pd(ai, b);       /* ai = [a.i*b.i a.i*b.r] */
    return _mm_addsub_pd(ar, ai); /* [a.r*b.r-a.i*b.i a.r*b.i+a.i*b.r] */
}

#define ZMUL(ar_, ai_, b_) \
    ar_ = ai_; \
    ar_ = _mm_movedup_pd(ar_); \
    ai_ = _mm_unpackhi_pd(ai_, ai_); \
    ar_ = _mm_mul_pd(ar_, b_); \
    b_  = _mm_shuffle_lh(b_, b_); \
    ai_ = _mm_mul_pd(ai_, b_); \
    ar_ = _mm_addsub_pd(ar_, ai_)


static void makeip(int n, int *ip);
static void makewt(int nw, int *ip, double *w);
static void makect(int nc, double *c);

static __inline void bitrv2 (int n, double *a, int const *ip);
static __inline void cftfsub(int n, double *a, double const *w);
static __inline void cftbsub(int n, double *a, double const *w);
static __inline void rftfsub(int n, double *a, double const *c);
static __inline void rftbsub(int n, double *a, double const *c);

static __inline void cft1st(int n, double *a, double const *w);
static __inline void cftmdl(int n, int l, double *a, double const *w);


void lsx_rdft_SSE3(int isgn, double *a, FFTcontext z)
{
    double xi;
    int n = z.n;

    if (isgn >= 0)
    {
        if (n > 4)
        {
            bitrv2(n, a, z.ip);
            cftfsub(n, a, z.wt);
            rftfsub(n, a, z.ct);
        }
        else if (n == 4)
        {
            cftfsub(n, a, z.wt);
        }
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    }
    else
    {
        a[1] = 0.5 * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4)
        {
            rftbsub(n, a, z.ct);
            bitrv2(n, a, z.ip);
            cftbsub(n, a, z.wt);
        }
        else if (n == 4)
        {
            cftfsub(n, a, z.wt);
        }
    }
}


/* -------- initializing routines -------- */

#include <stddef.h>
#include <malloc.h>

static __inline int dft_i_len(int len)
{
    int p = (int)(log((len) + .5) / M_LN2);
    int n = p/2 - 1;
    if (n < 0) n = 0;
    return 1 << n;
}

static __inline int dft_w_len(int len)
{
    return len;
}

static __inline int dft_c_len(int len)
{
    return len/2;
}


int lsx_rdft_init_SSE3(int n, FFTcontext* z)
{
    int nw = n >> 2;

    z->ip = NULL;
    z->wt = NULL;
    z->ct = NULL;
    z->n = n;

    if (n > 4)
    {
        z->ip = (int*  ) _aligned_malloc(dft_i_len(n) * sizeof(int),     16);
        if (!z->ip) { lsx_rdft_close(z); return -1; }
    }
    if (nw > 2)
    {
        z->wt = (double*) _aligned_malloc(dft_w_len(n) * sizeof(double), 16);
        if (!z->wt) { lsx_rdft_close(z); return -1; }
    }
    if (nw > 1)
    {
        z->ct = (double*) _aligned_malloc(dft_c_len(n) * sizeof(double), 16);
        if (!z->ct) { lsx_rdft_close(z); return -1; }
    }

    makewt(nw, z->ip, z->wt);
    makect(nw, z->ct);
    makeip(n, z->ip);

    return 0;
}


static void makeip(int n, int *ip)
{
    int j, l, m;

    if (n > 4)
        ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 3) < l) {
        l >>= 1;
        for (j = 0; j < m; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
}


static void makewt(int nw, int *ip, double *w)
{
    int j, nwh;
    double delta, x, y;

    if (nw > 2)
    {
        double* w0 = w + 3*nw;
        int k1, k2;
        nwh = nw >> 1;
        delta = M_PI_2 / nw;
        w0[0] = 1;
        w0[1] = 0;
        w0[nwh] = cos(M_PI_4);
        w0[nwh + 1] = w0[nwh];
        if (nwh > 2)
        {
            for (j = 2; j < nwh; j += 2)
            {
                x = cos(delta * j);
                y = sin(delta * j);
                w0[j  ] = x;
                w0[j+1] = y;
                w0[nw-j  ] = y;
                w0[nw-j+1] = x;
            }
            makeip(nw, ip);
            bitrv2(nw, w0, ip);
        }

        w[0] = w0[0]; w[1] = w0[1];
        k1 = 2; k2 = 4;
        for(j = 2; k2 < nw; j += 10)
        {
            w[j+0] = w0[k1+0];
            w[j+1] = w0[k1+1];

            w[j+2] = w0[k2+0];
            w[j+3] = w0[k2+1];

            w[j+4] = + w0[k2+0] - 2 * w0[k1+1] * w0[k2+1];
            w[j+5] = - w0[k2+1] + 2 * w0[k1+1] * w0[k2+0];

            w[j+6] = w0[k2+2];
            w[j+7] = w0[k2+3];

            w[j+8] = + w0[k2+2] - 2 * w0[k1] * w0[k2+3];
            w[j+9] = - w0[k2+3] + 2 * w0[k1] * w0[k2+2];

            k1 += 2; k2 += 4;
        }
    }
}


static void makect(int nc, double *c)
{
    int j, kk, nch, m;
    double delta, s1, c1;

    if (nc > 1) // nc >= 2, n >= 8, dft_c_len >= 4
    {
        nch = nc >> 1; //nch >= 1
        m = (nc << 1) + 1; // m >= 5
        delta = M_PI_2 / nc;

        c[0] = 0.0;
        c[1] = 0.0;
        c[nc  ] = 0.5 - 0.5 * M_SQRT1_2;
        c[nc+1] = 0.5 * M_SQRT1_2;

        kk = 2;
        for (j = 1; j < nch; j++)
        {
            s1 = 0.5 * sin(delta * j);
            c1 = 0.5 * cos(delta * j);

            c[kk    ] = 0.5 - s1;
            c[kk+1  ] = c1;
            c[m-kk-1] = 0.5 - c1;
            c[m-kk  ] = s1;

            kk += 2;
        }
    }
}


/* -------- child routines -------- */


static __inline void bitrv2(int n, double *a, int const *ip)
{
    int j, j1, k, k1, l, m, m2;
    //double xr, xi, yr, yi;

    l = n; m = 1;
    while ((m << 3) < l) { l >>= 1; m <<= 1; }

    m2 = 2 * m;
    if ((m << 3) == l) {
        for (k = 0; k < m; k++) {
            for (j = 0; j < k; j++) {
                __m128d x1, y1, x2, y2;
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x1 = _mm_load_pd(a+j1);
                y1 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y1);
                _mm_store_pd(a+k1, x1);

                j1 += m2;
                k1 += 2 * m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x2 = _mm_load_pd(a+j1);
                y2 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y2);
                _mm_store_pd(a+k1, x2);

                j1 += m2;
                k1 -= m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x1 = _mm_load_pd(a+j1);
                y1 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y1);
                _mm_store_pd(a+k1, x1);

                j1 += m2;
                k1 += 2 * m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x2 = _mm_load_pd(a+j1);
                y2 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y2);
                _mm_store_pd(a+k1, x2);
            }
            j1 = 2 * k + m2 + ip[k];
            k1 = j1 + m2;
            {
            __m128d x, y;
            //xr = a[j1];
            //xi = a[j1 + 1];
            //yr = a[k1];
            //yi = a[k1 + 1];
            x = _mm_load_pd(a+j1);
            y = _mm_load_pd(a+k1);
            //a[j1] = yr;
            //a[j1 + 1] = yi;
            //a[k1] = xr;
            //a[k1 + 1] = xi;
            _mm_store_pd(a+j1, y);
            _mm_store_pd(a+k1, x);
            }
        }
    }
    else
    {
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                __m128d x1, y1, x2, y2;
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x1 = _mm_load_pd(a+j1);
                y1 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y1);
                _mm_store_pd(a+k1, x1);

                j1 += m2;
                k1 += m2;

                //xr = a[j1];
                //xi = a[j1 + 1];
                //yr = a[k1];
                //yi = a[k1 + 1];
                x2 = _mm_load_pd(a+j1);
                y2 = _mm_load_pd(a+k1);
                //a[j1] = yr;
                //a[j1 + 1] = yi;
                //a[k1] = xr;
                //a[k1 + 1] = xi;
                _mm_store_pd(a+j1, y2);
                _mm_store_pd(a+k1, x2);
            }
        }
    }
}


static __inline void cftfsub(int n, double *a, double const *w)
{
    int j, j1, j2, j3, l;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        __m128d neg_hi = _mm_load_pd(C_PN);
        for (j = 0; j < l; j += 2) {
            __m128d x0, x1, x2, x3, t;
            j1 = j  + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r = a[j    ] + a[j1    ];
            //x0i = a[j + 1] + a[j1 + 1];
            //x1r = a[j    ] - a[j1    ];
            //x1i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            t  = _mm_load_pd(a+j1);
            x1 = x0;
            x1 = _mm_sub_pd(x1, t);
            x0 = _mm_add_pd(x0, t);

            //x2r = a[j2    ] + a[j3    ];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2    ] - a[j3    ];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t  = _mm_load_pd(a+j3);
            x3 = x2;
            x3 = _mm_sub_pd(x3, t);
            x2 = _mm_add_pd(x2, t);

            //a[j     ] = x0r + x2r;
            //a[j  + 1] = x0i + x2i;
            //a[j2    ] = x0r - x2r;
            //a[j2 + 1] = x0i - x2i;
            t  = x2;
            x2 = _mm_add_pd(x2, x0);
            x0 = _mm_sub_pd(x0, t);
            _mm_store_pd(a+j , x2);
            _mm_store_pd(a+j2, x0);

            //a[j1    ] = x1r - x3i;
            //a[j1 + 1] = x1i + x3r;
            //a[j3    ] = x1r + x3i;
            //a[j3 + 1] = x1i - x3r;
            x3 = _mm_shuffle_lh(x3, x3);
            x3 = _mm_neg_hi2(x3);
            t  = x3;
            x3 = _mm_add_pd(x3, x1);
            x1 = _mm_sub_pd(x1, t);
            _mm_store_pd(a+j1, x1);
            _mm_store_pd(a+j3, x3);
        }
    } else {
        for (j = 0; j < l; j += 2) {
            __m128d x0, x1, t;
            j1 = j + l;

            //x0r = a[j    ] - a[j1    ];
            //x0i = a[j + 1] - a[j1 + 1];
            //a[j    ] += a[j1    ];
            //a[j + 1] += a[j1 + 1];
            //a[j1    ] = x0r;
            //a[j1 + 1] = x0i;
            //    vv
            //a[j     ] = a[j    ] + a[j1    ];
            //a[j  + 1] = a[j + 1] + a[j1 + 1];
            //a[j1    ] = a[j    ] - a[j1    ];
            //a[j1 + 1] = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            x1 = _mm_load_pd(a+j1);
            t  = x1;
            x1 = _mm_add_pd(x1, x0);
            x0 = _mm_sub_pd(x0, t);
            _mm_store_pd(a+j , x1);
            _mm_store_pd(a+j1, x0);
        }
    }
}


static __inline void cftbsub(int n, double *a, double const *w)
{
    int j, j1, j2, j3, l;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        __m128d neg_hi = _mm_load_pd(C_PN);
        for (j = 0; j < l; j += 2) {
            __m128d x0, x1, x2, x3, t;
            j1 = j  + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r =  (a[j    ] + a[j1    ]);
            //x0i = -(a[j + 1] + a[j1 + 1]);
            //x1r =  (a[j    ] - a[j1    ]);
            //x1i = -(a[j + 1] - a[j1 + 1]);
            x0 = _mm_load_pd(a+j );
            t  = _mm_load_pd(a+j1);
            x1 = x0;
            x1 = _mm_sub_pd(x1, t);
            x0 = _mm_add_pd(x0, t);
            x1 = _mm_neg_hi2(x1);
            x0 = _mm_neg_hi2(x0);


            //x2r = a[j2    ] + a[j3    ];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2    ] - a[j3    ];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t  = _mm_load_pd(a+j3);
            x3 = x2;
            x3 = _mm_sub_pd(x3, t);
            x2 = _mm_add_pd(x2, t);

            //a[j     ] = x0r + x2r;
            //a[j  + 1] = x0i - x2i;
            //a[j2    ] = x0r - x2r;
            //a[j2 + 1] = x0i + x2i;
            x2 = _mm_neg_hi2(x2);
            t  = x0;
            x0 = _mm_sub_pd(x0, x2);
            x2 = _mm_add_pd(x2, t);
            _mm_store_pd(a+j , x2);
            _mm_store_pd(a+j2, x0);

            //a[j1    ] = x1r - x3i;
            //a[j1 + 1] = x1i - x3r;
            //a[j3    ] = x1r + x3i;
            //a[j3 + 1] = x1i + x3r;
            x3 = _mm_shuffle_lh(x3, x3);
            t  = x1;
            x1 = _mm_sub_pd(x1, x3);
            x3 = _mm_add_pd(x3, t);
            _mm_store_pd(a+j1, x1);
            _mm_store_pd(a+j3, x3);
        }
    } else {
        __m128d neg_hi = _mm_load_pd(C_PN);
        for (j = 0; j < l; j += 2) {
            __m128d x0, x1, t;
            j1 = j + l;

            //x0r =  (a[j    ] - a[j1    ]);
            //x0i = -(a[j + 1] - a[j1 + 1]);
            //...    vv
            //a[j     ] =  (a[j    ] + a[j1    ]);
            //a[j  + 1] = -(a[j + 1] + a[j1 + 1]);
            //a[j1    ] =  (a[j    ] - a[j1    ]);
            //a[j1 + 1] = -(a[j + 1] - a[j1 + 1]);
            x0 = _mm_load_pd(a+j );
            x1 = _mm_load_pd(a+j1);
            t  = x1;
            x1 = _mm_add_pd(x1, x0);
            x0 = _mm_sub_pd(x0, t);
            x1 = _mm_neg_hi2(x1);
            x0 = _mm_neg_hi2(x0);
            _mm_store_pd(a+j,  x1);
            _mm_store_pd(a+j1, x0);
        }
    }
}


static __inline void rftfsub(int n, double *a, double const *c)
{
    int j, k, kk, m;
    //double wkr, wki, xr, xi, yr, yi;
    __m128d neg_hi = _mm_load_pd(C_PN);

    //nc = n >> 2;
    m = n >> 1;
    kk = 2;
    for (j = 2; j < m; j += 2) {
        __m128d wk, aj, ak, x, y;
        k = n - j;

        //wkr = 0.5 - c[nc - kk];
        //wki = c[kk];
        wk = _mm_load_pd(c+kk);

        //xr = a[j    ] - a[k    ];
        //xi = a[j + 1] + a[k + 1];
        aj = _mm_load_pd(a+j);
        ak = _mm_load_pd(a+k);
        x = _mm_addsub_pd(aj, ak);

        //yr = wkr * xr - wki * xi;
        //yi = wkr * xi + wki * xr;
        ZMUL(y, wk, x);

        //a[j    ] = a[j    ] - yr;
        //a[j + 1] = a[j + 1] - yi;
        //a[k    ] = a[k    ] + yr;
        //a[k + 1] = a[k + 1] - yi;
        aj = _mm_sub_pd(aj, y);
        y  = _mm_neg_hi2(y);
        ak = _mm_add_pd(ak, y);
        _mm_store_pd(a+j, aj);
        _mm_store_pd(a+k, ak);

        kk += 2;
    }
}


static __inline void rftbsub(int n, double *a, double const *c)
{
    int j, k, kk, m;
    //double wkr, wki, xr, xi, yr, yi;
    __m128d neg_hi = _mm_load_pd(C_PN);

    //nc = n >> 2;
    m = n >> 1;
    kk = 2;
    a[1] = -a[1];
    for (j = 2; j < m; j += 2) {
        __m128d wk, aj, ak, x, y;
        k = n - j;

        //wkr = 0.5 - c[nc - kk];
        //wki = c[kk];
        wk = _mm_load_pd(c+kk);
        wk = _mm_neg_hi2(wk);

        //xr = a[j    ] - a[k    ];
        //xi = a[j + 1] + a[k + 1];
        aj = _mm_load_pd(a+j);
        ak = _mm_load_pd(a+k);
        x = _mm_addsub_pd(aj, ak);

        //yr = wkr * xr + wki * xi = wkr * xr - (-wki) * xi
        //yi = wkr * xi - wki * xr = wkr * xi + (-wki) * xr
        ZMUL(y, wk, x);

        //a[j    ] =  (a[j    ] - yr);
        //a[j + 1] = -(a[j + 1] - yi);
        //a[k    ] =  a[k    ] + yr;
        //a[k + 1] = -a[k + 1] + yi;
        aj = _mm_sub_pd(aj, y);
        aj = _mm_neg_hi2(aj);
        ak = _mm_neg_hi2(ak);
        ak = _mm_add_pd(ak, y);
        _mm_store_pd(a+j, aj);
        _mm_store_pd(a+k, ak);

        kk += 2;
    }
    a[m + 1] = -a[m + 1];
}


static __inline void cft1st(int n, double *a, double const *w)
{
    int j/*, k1, k2*/;
    //double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    __m128d x0, x1, x2, x3, t1, t2, wk2;
    __m128d neg_lo = _mm_load_pd(C_NP);

    //x0r = a[0] + a[2];
    //x0i = a[1] + a[3];
    //x1r = a[0] - a[2];
    //x1i = a[1] - a[3];
    x0 = _mm_load_pd(a+0);
    t1 = _mm_load_pd(a+2);
    x1 = x0;
    x0 = _mm_add_pd(x0, t1);
    x1 = _mm_sub_pd(x1, t1);

    //x2r = a[4] + a[6];
    //x2i = a[5] + a[7];
    //x3r = a[4] - a[6];
    //x3i = a[5] - a[7];
    x2 = _mm_load_pd(a+4);
    t2 = _mm_load_pd(a+6);
    x3 = x2;
    x2 = _mm_add_pd(x2, t2);
    x3 = _mm_sub_pd(x3, t2);

    //a[0] = x0r + x2r;
    //a[1] = x0i + x2i;
    //a[4] = x0r - x2r;
    //a[5] = x0i - x2i;
    t2 = x2;
    x2 = _mm_add_pd(x2, x0);
    x0 = _mm_sub_pd(x0, t2);
    _mm_store_pd(a  , x2);
    _mm_store_pd(a+4, x0);


    //a[2] = x1r - x3i;
    //a[3] = x1i + x3r;
    //a[6] = x1r + x3i;
    //a[7] = x1i - x3r;
    x3 = _mm_shuffle_lh(x3, x3);
    x3 = _mm_neg_lo2(x3);
    t1 = x3;
    x3 = _mm_add_pd(x3, x1);
    x1 = _mm_sub_pd(x1, t1);
    _mm_store_pd(a+2, x3);
    _mm_store_pd(a+6, x1);

    //x0r = a[8] + a[10];
    //x0i = a[9] + a[11];
    //x1r = a[8] - a[10];
    //x1i = a[9] - a[11];
    x0 = _mm_load_pd(a+8 );
    t1 = _mm_load_pd(a+10);
    x1 = x0;
    x0 = _mm_add_pd(x0, t1);
    x1 = _mm_sub_pd(x1, t1);

    //x2r = a[12] + a[14];
    //x2i = a[13] + a[15];
    //x3r = a[12] - a[14];
    //x3i = a[13] - a[15];
    x2 = _mm_load_pd(a+12);
    t2 = _mm_load_pd(a+14);
    x3 = x2;
    x2 = _mm_add_pd(x2, t2);
    x3 = _mm_sub_pd(x3, t2);

    //a[8] = x0r + x2r;
    //a[9] = x0i + x2i;
    //a[12] = x2i - x0i; = - (x0i - x2i)
    //a[13] = x0r - x2r; = + (x0r - x2r)
    t1 = x0;
    t1 = _mm_add_pd(t1, x2);
    x0 = _mm_sub_pd(x0, x2);
    x0 = _mm_shuffle_lh(x0, x0);
    x0 = _mm_neg_lo2(x0);
    _mm_store_pd(a+8 , t1);
    _mm_store_pd(a+12, x0);

    //x0r = x1r - x3i;
    //x0i = x1i + x3r;
    //x2r = x3i + x1r;
    //x2i = x3r - x1i;
    x3 = _mm_shuffle_lh(x3, x3);
    x0 = x1;
    x0 = _mm_addsub_pd(x0, x3);
    x1 = _mm_neg_lo2(x1);
    x3 = _mm_sub_pd(x3, x1);

    //wk1r = w[2];
    wk2 = _mm_load1_pd(w+2);

    //a[10] = wk1r * (x0r - x0i);
    //a[11] = wk1r * (x0r + x0i);
    //a[14] = wk1r * (x2i - x2r); = wk1r * (-x2r + x2i)
    //a[15] = wk1r * (x2i + x2r); = wk1r * (+x2r + x2i)
    t1 = x0;
    x0 = _mm_neg_hi(x0);
    x0 = _mm_hadd_pd(x0, t1);
    x2 = x3;
    x2 = _mm_neg_lo2(x2);
    x2 = _mm_hadd_pd(x2, x3);

    x0 = _mm_mul_pd(x0, wk2);
    x2 = _mm_mul_pd(x2, wk2);
    _mm_store_pd(a+10, x0);
    _mm_store_pd(a+14, x2);

    //k1 = 0; k2 = 0;
    w+=2;
    for (j = 16; j < n; j += 16)
    {
        //k1 += 2;
        //k2 += 4;
        //wk2r = w[k1];
        //wk2i = w[k1 + 1];
        //wk1r = w[k2];
        //wk1i = w[k2 + 1];
        //wk3r =   wk1r - 2 * wk2i * wk1i;
        //wk3i = - wk1i + 2 * wk2i * wk1r;

        wk2 = _mm_load_pd(w);
        //wk1 = _mm_load_pd(w+2);
        //wk3 = _mm_load_pd(w+4);

        //x0r = a[j    ] + a[j + 2];
        //x0i = a[j + 1] + a[j + 3];
        //x1r = a[j    ] - a[j + 2];
        //x1i = a[j + 1] - a[j + 3];
        x0 = _mm_load_pd(a+j+0);
        t1 = _mm_load_pd(a+j+2);
        x1 = x0;
        x1 = _mm_sub_pd(x1, t1);
        x0 = _mm_add_pd(x0, t1);

        //x2r = a[j + 4] + a[j + 6];
        //x2i = a[j + 5] + a[j + 7];
        //x3r = a[j + 4] - a[j + 6];
        //x3i = a[j + 5] - a[j + 7];
        x2 = _mm_load_pd(a+j+4);
        t2 = _mm_load_pd(a+j+6);
        x3 = x2;
        x3 = _mm_sub_pd(x3, t2);
        x2 = _mm_add_pd(x2, t2);

        //a[j    ] = x0r + x2r;
        //a[j + 1] = x0i + x2i;
        //x0r -= x2r;
        //x0i -= x2i;
        //a[j + 4] = wk2r * x0r - wk2i * x0i;
        //a[j + 5] = wk2r * x0i + wk2i * x0r;
        t1 = x2;
        x2 = _mm_add_pd(x2, x0);
        x0 = _mm_sub_pd(x0, t1);
        t1 = wk2;
        ZMUL(t2, t1, x0);
        _mm_store_pd(a+j  , x2);
        _mm_store_pd(a+j+4, t2);

        //x0r = x1r + x3i;
        //x0i = x1i - x3r;
        //x1r = x1r - x3i;
        //x1i = x1i + x3r;
        t1 = x3;
        t1 = _mm_neg_lo2(t1);
        t1 = _mm_shuffle_lh(t1, t1);
        x0 = x1;
        x0 = _mm_add_pd(x0, t1);
        x1 = _mm_sub_pd(x1, t1);

        //a[j + 2] = wk1r * x1r - wk1i * x1i;
        //a[j + 3] = wk1r * x1i + wk1i * x1r;
        //a[j + 6] = wk3r * x0r - wk3i * x0i;
        //a[j + 7] = wk3r * x0i + wk3i * x0r;
        x3 = _mm_load_pd(w+2);
        x2 = _mm_load_pd(w+4);
        ZMUL(t1, x3, x1);
        ZMUL(t2, x2, x0);
        _mm_store_pd(a+j+2, t1);
        _mm_store_pd(a+j+6, t2);

        //wk1r = w[k2 + 2];
        //wk1i = w[k2 + 3];
        //wk3r = + wk1r - 2 * wk2r * wk1i;
        //wk3i = - wk1i + 2 * wk2r * wk1r;
        //wk1 = _mm_load_pd(w+6);
        //wk3 = _mm_load_pd(w+8);

        //x0r = a[j + 8] + a[j + 10];
        //x0i = a[j + 9] + a[j + 11];
        //x1r = a[j + 8] - a[j + 10];
        //x1i = a[j + 9] - a[j + 11];
        x0 = _mm_load_pd(a+j+8);
        t1 = _mm_load_pd(a+j+10);
        x1 = x0;
        x1 = _mm_sub_pd(x1, t1);
        x0 = _mm_add_pd(x0, t1);

        //x2r = a[j + 12] + a[j + 14];
        //x2i = a[j + 13] + a[j + 15];
        //x3r = a[j + 12] - a[j + 14];
        //x3i = a[j + 13] - a[j + 15];
        x2 = _mm_load_pd(a+j+12);
        t2 = _mm_load_pd(a+j+14);
        x3 = x2;
        x3 = _mm_sub_pd(x3, t2);
        x2 = _mm_add_pd(x2, t2);

        //a[j + 8] = x0r + x2r;
        //a[j + 9] = x0i + x2i;
        //x0r -= x2r;
        //x0i -= x2i;
        //a[j + 12] = -wk2i * x0r - wk2r * x0i;
        //a[j + 13] = -wk2i * x0i + wk2r * x0r;
        t2 = x2;
        x2 = _mm_add_pd(x2, x0);
        x0 = _mm_sub_pd(x0, t2);
        ZMUL(t1, wk2, x0);
        t1 = _mm_shuffle_lh(t1, t1);
        t1 = _mm_neg_lo2(t1);
        _mm_store_pd(a+j+8, x2);
        _mm_store_pd(a+j+12, t1);

        //x0r = x1r + x3i;
        //x0i = x1i - x3r;
        //x1r = x1r - x3i;
        //x1i = x1i + x3r;
        t2 = x3;
        t2 = _mm_neg_lo2(t2);
        t2 = _mm_shuffle_lh(t2, t2);
        x0 = x1;
        x0 = _mm_add_pd(x0, t2);
        x1 = _mm_sub_pd(x1, t2);

        //a[j + 10] = wk1r * x1r - wk1i * x1i;
        //a[j + 11] = wk1r * x1i + wk1i * x1r;
        //a[j + 14] = wk3r * x0r - wk3i * x0i;
        //a[j + 15] = wk3r * x0i + wk3i * x0r;
        x3 = _mm_load_pd(w+6);
        x2 = _mm_load_pd(w+8);
        ZMUL(t1, x3, x1);
        ZMUL(t2, x2, x0);
        _mm_store_pd(a+j+10, t1);
        _mm_store_pd(a+j+14, t2);

        w+=10;
    }
}


static __inline void cftmdl(int n, int l, double *a, double const *w)
{
    int j, j1, j2, j3, k, /*k1, k2,*/ m, m2;
    //double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    //double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    __m128d wk1;

    m = l << 2;
    for (j = 0; j < l; j += 2)
    {
        __m128d x0, x1, x2, x3, t1, t2;
        __m128d neg_lo = _mm_load_pd(C_NP);
        j1 = j  + l;
        j2 = j1 + l;
        j3 = j2 + l;

        //x0r = a[j    ] + a[j1    ];
        //x0i = a[j + 1] + a[j1 + 1];
        //x1r = a[j    ] - a[j1    ];
        //x1i = a[j + 1] - a[j1 + 1];
        x0 = _mm_load_pd(a+j );
        t1 = _mm_load_pd(a+j1);
        x1 = x0;
        x1 = _mm_sub_pd(x1, t1);
        x0 = _mm_add_pd(x0, t1);

        //x2r = a[j2    ] + a[j3    ];
        //x2i = a[j2 + 1] + a[j3 + 1];
        //x3r = a[j2    ] - a[j3    ];
        //x3i = a[j2 + 1] - a[j3 + 1];
        x2 = _mm_load_pd(a+j2);
        t2 = _mm_load_pd(a+j3);
        x3 = x2;
        x3 = _mm_sub_pd(x3, t2);
        x2 = _mm_add_pd(x2, t2);

        //a[j    ] = x0r + x2r;
        //a[j + 1] = x0i + x2i;
        //a[j2    ] = x0r - x2r;
        //a[j2 + 1] = x0i - x2i;
        t1 = x2;
        x2 = _mm_add_pd(x2, x0);
        x0 = _mm_sub_pd(x0, t1);
        _mm_store_pd(a+j , x2);
        _mm_store_pd(a+j2, x0);

        //a[j1    ] = x1r - x3i;
        //a[j1 + 1] = x1i + x3r;
        //a[j3    ] = x1r + x3i;
        //a[j3 + 1] = x1i - x3r;
        x3 = _mm_shuffle_lh(x3, x3);
        x3 = _mm_neg_lo2(x3);
        t2 = x3;
        x3 = _mm_add_pd(x3, x1);
        x1 = _mm_sub_pd(x1, t2);
        _mm_store_pd(a+j1, x3);
        _mm_store_pd(a+j3, x1);
    }
    //wk1r = w[2];
    wk1 = _mm_load1_pd(w+2);
    {
        __m128d neg_lo = _mm_load_pd(C_NP), neg_hi = _mm_load_pd(C_PN);
        for (j = m; j < l + m; j += 2)
        {
            __m128d x0, x1, x2, x3, t;
            j1 = j  + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r = a[j] + a[j1];
            //x0i = a[j + 1] + a[j1 + 1];
            //x1r = a[j] - a[j1];
            //x1i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            t  = _mm_load_pd(a+j1);
            x1 = x0;
            x1 = _mm_sub_pd(x1, t);
            x0 = _mm_add_pd(x0, t);

            //x2r = a[j2] + a[j3];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2] - a[j3];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t  = _mm_load_pd(a+j3);
            x3 = x2;
            x3 = _mm_sub_pd(x3, t);
            x2 = _mm_add_pd(x2, t);

            //a[j    ] = x0r + x2r;
            //a[j + 1] = x0i + x2i;
            //a[j2    ] = x2i - x0i; = - (x0i - x2i)
            //a[j2 + 1] = x0r - x2r; = + (x0r - x2r)

            t  = x2;
            x2 = _mm_add_pd(x2, x0);
            x0 = _mm_sub_pd(x0, t);
            x0 = _mm_shuffle_lh(x0, x0);
            x0 = _mm_neg_lo2(x0);
            _mm_store_pd(a+j , x2);
            _mm_store_pd(a+j2, x0);

            //x0r = x1r - x3i;
            //x0i = x1i + x3r;
            //a[j1    ] = wk1r * (x0r - x0i);
            //a[j1 + 1] = wk1r * (x0r + x0i);
            //x0r = x3i + x1r;
            //x0i = x3r - x1i;
            //a[j3    ] = wk1r * (x0i - x0r); = wk1r * (- x0r + x0i)
            //a[j3 + 1] = wk1r * (x0i + x0r); = wk1r * (+ x0r + x0i)
            x3 = _mm_shuffle_lh(x3, x3);

            x0 = x1;
            x0 = _mm_addsub_pd(x0, x3);
            x2 = x0;
            x0 = _mm_neg_hi2(x0);
            x0 = _mm_hadd_pd(x0, x2);

            x1 = _mm_neg_lo2(x1);
            x3 = _mm_sub_pd(x3, x1);
            x2 = x3;
            x3 = _mm_neg_lo2(x3);
            x3 = _mm_hadd_pd(x3, x2);

            x0 = _mm_mul_pd(x0, wk1);
            x3 = _mm_mul_pd(x3, wk1);
            _mm_store_pd(a+j1, x0);
            _mm_store_pd(a+j3, x3);
        }
    }

    //k1 = 0; k2 = 0;
    m2 = 2 * m;
    w+=2;
    for (k = m2; k < n; k += m2)
    {
        __m128d x0, x1, x2, x3, wk2, wk3, t;
        //k1 += 2;
        //k2 += 4;

        //wk2r = w[k1];
        //wk2i = w[k1 + 1];
        //wk1r = w[k2];
        //wk1i = w[k2 + 1];
        //wk3r = wk1r - 2 * wk2i * wk1i; = + (wk1r - 2 * wk2i * wk1i)
        //wk3i = 2 * wk2i * wk1r - wk1i; = - (wk1i - 2 * wk2i * wk1r)

        wk2 = _mm_load_pd(w  );
        wk1 = _mm_load_pd(w+2);
        wk3 = _mm_load_pd(w+4);

        for (j = k; j < l + k; j += 2)
        {
            j1 = j  + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r = a[j    ] + a[j1    ];
            //x0i = a[j + 1] + a[j1 + 1];
            //x1r = a[j    ] - a[j1    ];
            //x1i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            t  = _mm_load_pd(a+j1);
            x1 = x0;
            x1 = _mm_sub_pd(x1, t);
            x0 = _mm_add_pd(x0, t);
            //x2r = a[j2] + a[j3];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2] - a[j3];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t  = _mm_load_pd(a+j3);
            x3 = x2;
            x3 = _mm_sub_pd(x3, t);
            x2 = _mm_add_pd(x2, t);

            //a[j    ] = x0r + x2r;
            //a[j + 1] = x0i + x2i;
            //x0r -= x2r;
            //x0i -= x2i;
            //a[j2    ] = wk2r * x0r - wk2i * x0i;
            //a[j2 + 1] = wk2r * x0i + wk2i * x0r;
            t  = x2;
            x2 = _mm_add_pd(x2, x0);
            x0 = _mm_sub_pd(x0, t);
            t  = zmul(wk2, x0);
            _mm_store_pd(a+j, x2);
            _mm_store_pd(a+j2, t);

            //x0r = x1r + x3i;
            //x0i = x1i - x3r;
            //x1r = x1r - x3i;
            //x1i = x1i + x3r;
            t = _mm_neg_lo(x3);
            t = _mm_shuffle_lh(t, t);
            x0 = x1;
            x0 = _mm_add_pd(x0, t);
            x1 = _mm_sub_pd(x1, t);

            //a[j1    ] = wk1r * x1r - wk1i * x1i;
            //a[j1 + 1] = wk1r * x1i + wk1i * x1r;
            //a[j3    ] = wk3r * x0r - wk3i * x0i;
            //a[j3 + 1] = wk3r * x0i + wk3i * x0r;
            _mm_store_pd(a+j3, zmul(wk3, x0));
            _mm_store_pd(a+j1, zmul(wk1, x1));
        }

        //wk1r = w[k2 + 2];
        //wk1i = w[k2 + 3];
        //wk3r = + wk1r - 2 * wk2r * wk1i;
        //wk3i = - wk1i + 2 * wk2r * wk1r;

        wk1 = _mm_load_pd(w+6);
        wk3 = _mm_load_pd(w+8);

        for (j = k + m; j < l + (k + m); j += 2)
        {
            j1 = j  + l;
            j2 = j1 + l;
            j3 = j2 + l;

            //x0r = a[j    ] + a[j1    ];
            //x0i = a[j + 1] + a[j1 + 1];
            //x1r = a[j    ] - a[j1    ];
            //x1i = a[j + 1] - a[j1 + 1];
            x0 = _mm_load_pd(a+j );
            t  = _mm_load_pd(a+j1);
            x1 = x0;
            x1 = _mm_sub_pd(x1, t);
            x0 = _mm_add_pd(x0, t);

            //x2r = a[j2    ] + a[j3    ];
            //x2i = a[j2 + 1] + a[j3 + 1];
            //x3r = a[j2    ] - a[j3    ];
            //x3i = a[j2 + 1] - a[j3 + 1];
            x2 = _mm_load_pd(a+j2);
            t  = _mm_load_pd(a+j3);
            x3 = x2;
            x3 = _mm_sub_pd(x3, t);
            x2 = _mm_add_pd(x2, t);

            //a[j    ] = x0r + x2r;
            //a[j + 1] = x0i + x2i;
            //x0r -= x2r;
            //x0i -= x2i;
            //a[j2    ] = -wk2i * x0r - wk2r * x0i;
            //a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
            t = x2;
            x2 = _mm_add_pd(x2, x0);
            x0 = _mm_sub_pd(x0, t);
            t = zmul(wk2, x0);
            t = _mm_shuffle_lh(t, t);
            t = _mm_neg_lo(t);
            _mm_store_pd(a+j, x2);
            _mm_store_pd(a+j2, t);

            //x0r = x1r + x3i;
            //x0i = x1i - x3r;
            //x1r = x1r - x3i;
            //x1i = x1i + x3r;
            t = _mm_neg_lo(x3);
            t = _mm_shuffle_lh(t, t);
            x0 = x1;
            x0 = _mm_add_pd(x0, t);
            x1 = _mm_sub_pd(x1, t);
            //a[j1    ] = wk1r * x1r - wk1i * x1i;
            //a[j1 + 1] = wk1r * x1i + wk1i * x1r;
            //a[j3    ] = wk3r * x0r - wk3i * x0i;
            //a[j3 + 1] = wk3r * x0i + wk3i * x0r;
            _mm_store_pd(a+j1, zmul(wk1, x1));
            _mm_store_pd(a+j3, zmul(wk3, x0));
        }
        w += 10;
    }
}
