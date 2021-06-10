/* Internal interface to SoX rate() resampling routines
 *
 * Copyright (c) 2011-2012 lvqcl.
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

#ifndef RATE_I_H
#define RATE_I_H

#include "ratelib.h"

typedef struct RR_vtable_tag
{
    int  (*init)(RR_handle* h, const RR_config* config, int nchannels); /*create+start*/
    int  (*flow)(RR_handle *h, const fb_sample_t* ibuf, fb_sample_t* obuf, size_t isamp, size_t osamp, size_t* iused, size_t* ogen);
    int  (*push)(RR_handle *h, const fb_sample_t* ibuf, size_t isamp);
    int  (*pull)(RR_handle *h, fb_sample_t* obuf, size_t osamp, size_t* ogen);
    int  (*drain)(RR_handle *h);
    void (*close)(RR_handle *h); /*stop*/
} RR_vtable;

#ifdef __cplusplus
extern "C" {
#endif

RR_handle* RR_ctor_SSE3  (const RR_config* config, int nchannels);
RR_handle* RR_ctor_double(const RR_config* config, int nchannels);

RR_handle* RR_ctor_SSE   (const RR_config* config, int nchannels);
RR_handle* RR_ctor_float (const RR_config* config, int nchannels);

#ifdef __cplusplus
}
#endif

#endif
