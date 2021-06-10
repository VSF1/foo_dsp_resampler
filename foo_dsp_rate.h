/* Copyright (c) 2008-12, 2018 lvqcl. All rights reserved.
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

#ifndef FOODSPRATE_H
#define FOODSPRATE_H

#include "dsp_config.h"
#include "chain.h"

class dsp_rate : public dsp_impl_base_t<dsp_v2>
{
    resampler_link rate_;
    RateConfig cfg_;

    size_t in_samples_accum_, out_samples_gen_accum_; // for latency

    audio_sample* in_buffer_0;
    audio_sample* in_buffer_;
    audio_sample* out_buffer_;
    unsigned BUF_CHANNELS_;
    size_t INBUF_SIZE_0;
    size_t INBUF_SIZE_;
    size_t OUTBUF_SIZE_;
    size_t PRIME_LEN_;

    unsigned out_rate_;
    unsigned sample_rate_;
    unsigned channel_count_;
    unsigned channel_map_;

    size_t N_samples_to_add_, N_samples_to_drop_;
    size_t samples_in_buffer_;
    size_t samples_dropped_;
    bool is_preextrapolated_;

    void ctor_init();
    bool set_data(const dsp_preset & p_data);

public:
    dsp_rate();
    dsp_rate(const t_dsp_rate_params& params);
    ~dsp_rate();

private:
    void reinit(unsigned sample_rate, unsigned channel_count, unsigned channel_map);
    void close();
    void flushwrite();

    void as_memcpy2 (audio_sample* dst0, size_t off_d, audio_sample* src0, size_t off_s, size_t samples)
        { memcpy (dst0 + off_d*channel_count_, src0 + off_s*channel_count_, samples * channel_count_ * sizeof(audio_sample)); }

    void as_memmove2(audio_sample* dst0, size_t off_d, audio_sample* src0, size_t off_s, size_t samples)
        { memmove(dst0 + off_d*channel_count_, src0 + off_s*channel_count_, samples * channel_count_ * sizeof(audio_sample)); }

    static const size_t LPC_ORDER = 32;

protected:
    virtual void on_endoftrack(abort_callback & p_abort);
    virtual void on_endofplayback(abort_callback & p_abort);
    virtual bool on_chunk(audio_chunk * chunk, abort_callback & p_abort);

public:
    virtual void flush();
    virtual double get_latency();
    virtual bool need_track_change_mark() { return false; }
};

#endif
