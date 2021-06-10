/* Interface to SoX rate() resampling routines
 *
 * Copyright (c) 2008-11 lvqcl.
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

#ifndef RATELIB_H
#define RATELIB_H

#include <stddef.h>

enum RR_error
{
    RR_OK = 0,
    RR_ENOMEM,
    RR_INTERNAL,
    RR_NULLHANDLE,
    RR_RATEERROR,
    RR_EXTUNINIT,
    RR_INVPARAM,
};

enum RR_quality
{
    //RR_default = -1,

    RR_best = 0,
    RR_norm = 1,
};

enum RR_phase
{
    RR_minimum = 0,
    RR_linear  = 50,
    RR_maximum = 100,
};

typedef float fb_sample_t;

typedef struct RR_config_tag
{
    size_t in_rate;
    size_t out_rate;

    double phase;
    double bandwidth;
    int allow_aliasing;

    enum RR_quality quality;
} RR_config;

typedef struct RR_handle_tag RR_handle;


#ifdef __cplusplus
extern "C" {
#endif

int init_ratelib(void (*alloc_error_handler)(void));

int RR_open(const RR_config *config, int nchannels, RR_handle **const handle);
int RR_flow(RR_handle *h, const fb_sample_t *ibuf, fb_sample_t *obuf, size_t isamp, size_t osamp, size_t *iused, size_t *ogen);
int RR_push(RR_handle *h, const fb_sample_t *ibuf, size_t isamp);
int RR_pull(RR_handle *h, fb_sample_t *obuf, size_t osamp, size_t *ogen);
int RR_drain(RR_handle *h);
void RR_close(RR_handle **h);

const char* RR_strerror(int error);

#ifdef __cplusplus
}
#endif

#endif
