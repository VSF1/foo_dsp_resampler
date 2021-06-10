/* Copyright (c) 2018 lvqcl.  All rights reserved.
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

#ifndef CHAIN_H_
#define CHAIN_H_

#include "rate/ratelib.h"
class resampler_link
{
    RR_handle * handle;

    static void msg(int err)
    {
        if (err) { console::error(RR_strerror(err)); throw std::bad_alloc(); }
    }

public:
    resampler_link() : handle(NULL) {}
    ~resampler_link() { if (handle) close(); }
    bool is_initialized() { return (handle != NULL); }

    void open(const RR_config * cfg, unsigned nch) { int err = RR_open(cfg, nch, &handle); msg(err); }
    void push(const audio_sample *ibuf, size_t isamp) { int err = RR_push(handle, ibuf, isamp); msg(err); }
    void pull(audio_sample *obuf, size_t osamp, size_t *ogen) { int err = RR_pull(handle, obuf, osamp, ogen); msg(err); }
    void drain() { int err = RR_drain(handle); msg(err); }
    void close() { RR_close(&handle); } // zeroes handle

    /* void flow(const audio_sample *ibuf, audio_sample *obuf, size_t isamp, size_t osamp, size_t *iused, size_t *ogen) { int err = RR_flow(handle, ibuf, obuf, isamp, osamp, iused, ogen); msg(err); } */
};

#if 0
class copy_link
{
    audio_sample* buf;
    size_t bufsize;
    size_t in_buf;
    unsigned nch;

    void cc_memcpy(audio_sample* dst, const audio_sample* src, size_t samples) { memcpy(dst, src, samples * nch * sizeof(audio_sample)); }
    void cc_memmove(audio_sample* dst, const audio_sample* src, size_t samples) { memmove(dst, src, samples * nch * sizeof(audio_sample)); }

public:
    copy_link() : buf(NULL), bufsize(0), in_buf(0), nch(0) {}
    ~copy_link() { if (buf) close(); }
    bool is_initialized() { return (buf != NULL); }

    void open(const RR_config * cfg, unsigned nchan) { nch = nchan; }
    void drain() {}
    void close() { delete[] buf; buf = 0; }

    void push(const audio_sample *ibuf, size_t isamp)
    {
        size_t avail = bufsize - in_buf;
        if (avail < isamp)
        {
            size_t newsize = in_buf + isamp + isamp/2;
            audio_sample* newbuf = new audio_sample[newsize*nch];
            cc_memcpy(newbuf, buf, in_buf);
            delete[] buf;
            buf = newbuf; bufsize = newsize;
        }
        cc_memcpy(buf + in_buf*nch, ibuf, isamp);
        in_buf += isamp;
    }

    void pull(audio_sample *obuf, size_t osamp, size_t *ogen)
    {
        size_t to_send = min(in_buf, osamp);
        if (ogen) *ogen = to_send;

        cc_memcpy(obuf, buf, to_send);
        in_buf -= to_send;
        cc_memmove(buf, buf + to_send*nch, in_buf);
    }
};
#endif

#endif
