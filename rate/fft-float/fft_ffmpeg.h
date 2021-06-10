/*
 * Copyright (c) 2000, 2001, 2002 Fabrice Bellard
 * Copyright (c) 2002-2004 Michael Niedermayer <michaelni@gmx.at>
 * Copyright (c) 2009 Alex Converse <alex dot converse at gmail dot com>
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef FFT_H
#define FFT_H

#include <stdint.h>
#define FF_ALIGN 32

#ifdef __cplusplus
extern "C" {
#endif

static __inline void *ff_malloc(unsigned int size)
{
    return _aligned_malloc(size, FF_ALIGN);
}

static __inline void ff_free(void *ptr)
{
    _aligned_free(ptr);
}

static __inline void ff_freep(void* *ptr)
{
    ff_free(*ptr);
    *ptr = NULL;
}



typedef float FFTSample;

typedef struct FFTComplex {
    FFTSample re, im;
} FFTComplex;

typedef struct FFTContext FFTContext;

struct FFTContext {
    int nbits;
    int inverse;
    uint16_t *revtab;
};

void ff_fft_init_costables();

int ff_fft_init(FFTContext *s, int nbits, int inverse, int sse);
void ff_fft_end(FFTContext *s);
/**
 * Do the permutation needed BEFORE calling ff_fft_calc().
 */
void ff_fft_permute_c(FFTContext *s, FFTComplex *z, FFTComplex *tmp_buf);
void ff_fft_permute_sse(FFTContext *s, FFTComplex *z, FFTComplex *tmp_buf);
/**
 * Do a complex FFT with the parameters defined in ff_fft_init(). The
 * input data must be permuted before. No 1.0/sqrt(n) normalization is done.
 */
void ff_fft_calc_c(FFTContext *s, FFTComplex *z);
void ff_fft_calc_sse(FFTContext *s, FFTComplex *z);


enum RDFTransformType {
    DFT_R2C,
    IDFT_C2R,
    IDFT_R2C,
    DFT_C2R,
};

typedef struct RDFTContext RDFTContext;

struct RDFTContext {
    int nbits;
    int inverse;
    int sign_convention;

    /* pre/post rotation tables */
    const FFTSample *tcos;
    const FFTSample *tsin;
    int negative_sin;
    FFTContext fft;
};

int ff_rdft_init(RDFTContext *s, int nbits, enum RDFTransformType trans, int sse);
void ff_rdft_end(RDFTContext *s);

void ff_rdft_calc_c(RDFTContext *s, FFTSample *data, FFTComplex *tmp_buf);
void ff_rdft_calc_sse(RDFTContext *s, FFTSample *data, FFTComplex *tmp_buf);

#ifdef __cplusplus
}
#endif

#endif
