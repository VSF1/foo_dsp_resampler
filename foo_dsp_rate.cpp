/* Copyright (c) 2008-12, 18 lvqcl. All rights reserved.
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

#include "stdafx.h"
#include "new_exception.h"
#include "lpc/lpc.h"
#include "util.h"

__declspec(noreturn) static inline void throw_new_exception_() { throw std::bad_alloc(); }
extern "C" void throw_new_exception() { throw_new_exception_(); }


void dsp_rate::ctor_init()
{
    t_dsp_rate_params params;
    params.get_rateconfig(cfg_);

    in_samples_accum_ = out_samples_gen_accum_ = 0;

    in_buffer_0 = NULL; in_buffer_ = NULL;
    INBUF_SIZE_0 = 0; INBUF_SIZE_ = 0;
    out_buffer_ = NULL;
    OUTBUF_SIZE_ = 0;
    BUF_CHANNELS_ = 0;
    PRIME_LEN_ = 0;

    out_rate_ = 0;
    sample_rate_ = 0;
    channel_count_ = 0;
    channel_map_ = 0;

    N_samples_to_add_ = N_samples_to_drop_ = 0;

    is_preextrapolated_ = false;
    samples_in_buffer_ = 0;
    samples_dropped_ = 0;
}

dsp_rate::dsp_rate()
{
    ctor_init();
}

dsp_rate::dsp_rate(const t_dsp_rate_params & params)
{
    ctor_init();
    params.get_rateconfig(cfg_);
}

bool dsp_rate::set_data(const dsp_preset & p_data)
{
    t_dsp_rate_params params;
    if (!params.set_data(p_data)) return false;
    params.get_rateconfig(cfg_);
    return true;
}

dsp_rate::~dsp_rate()
{
    close(); // rate_.is_initialised() == false

    delete[] in_buffer_0; in_buffer_0 = NULL; INBUF_SIZE_0 = 0;
    delete[] out_buffer_; out_buffer_ = NULL; OUTBUF_SIZE_ = 0; BUF_CHANNELS_ = 0;
}

void dsp_rate::on_endoftrack(abort_callback & p_abort) { flushwrite(); }

void dsp_rate::on_endofplayback(abort_callback & p_abort) { flushwrite(); }

void dsp_rate::reinit(unsigned sample_rate, unsigned channel_count, unsigned channel_map)
{
    out_rate_ = cfg_.realrate(sample_rate);
    RR_config c = { sample_rate, out_rate_, double(cfg_.phase), double(cfg_.passband10)/10.0, cfg_.allow_aliasing?1:0, (RR_quality)cfg_.quality };

    rate_.open(&c, channel_count);

    channel_count_ = channel_count; channel_map_ = channel_map; sample_rate_ = sample_rate;
    in_samples_accum_ = out_samples_gen_accum_ = 0;

    samples_in_buffer_ = 0; samples_dropped_ = 0;
    is_preextrapolated_ = false;
    N_samples_to_add_ = sample_rate;
    N_samples_to_drop_ = out_rate_;

    samples_len(&N_samples_to_add_, &N_samples_to_drop_, 20, 8192u);
    INBUF_SIZE_ = max(sample_rate/10, 2048u); INBUF_SIZE_ = min(INBUF_SIZE_, 65536u);
    PRIME_LEN_ = max(sample_rate/20, 1024u); PRIME_LEN_  = min(PRIME_LEN_ , 16384u); PRIME_LEN_ = max(PRIME_LEN_, 2*LPC_ORDER + 1);

    size_t new_inbuf_size = max(INBUF_SIZE_0, N_samples_to_add_ + INBUF_SIZE_ + N_samples_to_add_);
    size_t new_outbuf_size = max(OUTBUF_SIZE_, min(next_pow2(out_rate_/10, 8192u), 65536u) + N_samples_to_drop_);
    unsigned new_buf_channels = max(BUF_CHANNELS_, channel_count);

    if (BUF_CHANNELS_ < new_buf_channels || INBUF_SIZE_0 < new_inbuf_size || OUTBUF_SIZE_ < new_outbuf_size)
    {
        delete[] in_buffer_0;
        delete[] out_buffer_;

        in_buffer_0 = new audio_sample[new_buf_channels*new_inbuf_size];
        out_buffer_ = new audio_sample[new_buf_channels*new_outbuf_size];

        INBUF_SIZE_0 = new_inbuf_size;
        OUTBUF_SIZE_ = new_outbuf_size;
        BUF_CHANNELS_ = new_buf_channels;
    }

    in_buffer_ = in_buffer_0 + N_samples_to_add_*channel_count_;
}

void dsp_rate::close()
{
    rate_.close(); // rate_.is_initialised() == false

    in_samples_accum_ = out_samples_gen_accum_ = 0; // necessary?
}

bool dsp_rate::on_chunk(audio_chunk * chunk, abort_callback & p_abort)
{
    unsigned channel_count = chunk->get_channels();
    unsigned channel_map = chunk->get_channel_config();
    unsigned sample_rate = chunk->get_sample_rate();
    t_size sample_count = chunk->get_sample_count();
    audio_sample * current = chunk->get_data();

    if (!rate_.is_initialized())
    {
        if (cfg_.is_no_resample(sample_rate)) return true;
        reinit(sample_rate, channel_count, channel_map);
    }
    else if ((channel_count_ != channel_count) || (channel_map_ != channel_map) || (sample_rate_ != sample_rate))
    {    // number of channels or samplerate has changed - reinitialize
        flushwrite(); //here channel_count_, channel_map_ and sample_rate_ must have old values
        // link.is_initialized() == false here
        if (cfg_.is_no_resample(sample_rate)) return true;
        reinit(sample_rate, channel_count, channel_map);
    }

    size_t out_samples_gen;
    do {

        if (!is_preextrapolated_)
        {
            size_t samples_to_move = min(sample_count, INBUF_SIZE_ - samples_in_buffer_);
            as_memcpy2(in_buffer_, samples_in_buffer_, current, 0, samples_to_move);
            samples_in_buffer_ += samples_to_move;
            current += samples_to_move*channel_count_;
            sample_count -=samples_to_move;
            in_samples_accum_ += samples_to_move;

            if (samples_in_buffer_ == INBUF_SIZE_)
            {
                lpc_extrapolate_bkwd(in_buffer_, INBUF_SIZE_, PRIME_LEN_, channel_count_, LPC_ORDER, N_samples_to_add_);
                is_preextrapolated_ = true;
                rate_.push(in_buffer_0, N_samples_to_add_ + INBUF_SIZE_);
            }
        }

        if (is_preextrapolated_ && sample_count)    // samples_in_buffer_ == INBUF_SIZE_
        {
            if (sample_count < INBUF_SIZE_)
            {
                as_memmove2(in_buffer_, 0, in_buffer_, sample_count, INBUF_SIZE_ - sample_count);
                as_memcpy2(in_buffer_, INBUF_SIZE_ - sample_count, current, 0, sample_count);
            }
            else // sample_count >= INBUF_SIZE_
            {
                as_memcpy2(in_buffer_, 0, current, sample_count - INBUF_SIZE_, INBUF_SIZE_);
            }
            rate_.push(current, sample_count);
            current += sample_count*channel_count_;
            in_samples_accum_ += sample_count; sample_count = 0;
        }

        rate_.pull(out_buffer_, OUTBUF_SIZE_, &out_samples_gen);

        size_t to_drop = N_samples_to_drop_ - samples_dropped_;
        if (to_drop)
        {
            to_drop = min(to_drop, out_samples_gen);
            out_samples_gen -= to_drop;
            samples_dropped_ += to_drop;
        }
        if (out_samples_gen)
        {
            out_samples_gen_accum_ += out_samples_gen;
            audio_chunk * out = insert_chunk(out_samples_gen*channel_count_);
            out->set_data(out_buffer_ + to_drop*channel_count_, out_samples_gen, channel_count_, out_rate_, channel_map_);
        }
    } while (sample_count || out_samples_gen);

    while (in_samples_accum_ > sample_rate_ && out_samples_gen_accum_ > out_rate_)
    {
        in_samples_accum_ -= sample_rate_;
        out_samples_gen_accum_ -= out_rate_;
    }
    return false;
}

void dsp_rate::flush()
{
    if (!rate_.is_initialized()) return;
    close();
}

void dsp_rate::flushwrite()
{
    if (!rate_.is_initialized()) return;

    if (!is_preextrapolated_ && !(samples_in_buffer_ > 2 * LPC_ORDER))    // cannot extrapolate
    {
        rate_.push(in_buffer_, samples_in_buffer_);
        rate_.drain();

        size_t out_samples_gen;
        while(1)
        {
            rate_.pull(out_buffer_, OUTBUF_SIZE_, &out_samples_gen);
            if (out_samples_gen == 0) break;

            out_samples_gen_accum_ += out_samples_gen;
            audio_chunk * out = insert_chunk(out_samples_gen*channel_count_);
            out->set_data(out_buffer_, out_samples_gen, channel_count_, out_rate_, channel_map_);
        }
        close();
        return;
    }

    if (!is_preextrapolated_)    // samples_in_buffer_ > 2 * LPC_ORDER, can extrapolate
    {
        size_t prime = min(samples_in_buffer_, PRIME_LEN_);
        lpc_extrapolate_bkwd(in_buffer_, samples_in_buffer_, prime, channel_count_, LPC_ORDER, N_samples_to_add_);
        lpc_extrapolate_fwd (in_buffer_, samples_in_buffer_, prime, channel_count_, LPC_ORDER, N_samples_to_add_);
        is_preextrapolated_ = true;
        rate_.push(in_buffer_0, N_samples_to_add_ + samples_in_buffer_ + N_samples_to_add_);
        rate_.drain();

        samples_dropped_ = 0;
        samples_in_buffer_ = 0;
        size_t out_samples_gen;

        while (1)
        {
            rate_.pull(out_buffer_ + samples_in_buffer_*channel_count_, OUTBUF_SIZE_ - samples_in_buffer_, &out_samples_gen);
            if (out_samples_gen == 0) break;

            // now to drop first (N_samples_to_drop_) samples...
            size_t to_drop = N_samples_to_drop_ - samples_dropped_;
            to_drop = min(to_drop, out_samples_gen);
            if (to_drop)
            {
                // assert(samples_in_buffer_ == 0);
                out_samples_gen -= to_drop;
                samples_dropped_ += to_drop;
                as_memmove2(out_buffer_, 0, out_buffer_, to_drop, out_samples_gen);
            }

            // ...and also last (N_samples_to_drop_) samples;
            samples_in_buffer_ += out_samples_gen;
            size_t avail = samples_in_buffer_ - min(samples_in_buffer_, N_samples_to_drop_);
            if (avail)
            {
                out_samples_gen_accum_ += avail;
                audio_chunk * out = insert_chunk(avail*channel_count_);
                out->set_data(out_buffer_, avail, channel_count_, out_rate_, channel_map_);
                samples_in_buffer_ -= avail; /* samples_in_buffer_ = N_samples_to_drop_ */
                as_memmove2(out_buffer_, 0, out_buffer_, avail, samples_in_buffer_);
            }
        }
        close();
        return;
    }

    // is_preextrapolated_ == true; need to post-extrapolate; samples_in_buffer_ == INBUF_SIZE_, but they are already processed
    {
        lpc_extrapolate_fwd(in_buffer_, INBUF_SIZE_, PRIME_LEN_, channel_count_, LPC_ORDER, N_samples_to_add_);
        rate_.push(in_buffer_ + INBUF_SIZE_*channel_count_, N_samples_to_add_); // process N added samples
        rate_.drain();

        samples_in_buffer_ = 0;
        size_t out_samples_gen;
        while (1)
        {
            rate_.pull(out_buffer_ + samples_in_buffer_*channel_count_, OUTBUF_SIZE_ - samples_in_buffer_, &out_samples_gen);
            if (out_samples_gen == 0) break; /* assert (samples_in_buffer_ == N_samples_to_drop_) */

            samples_in_buffer_ += out_samples_gen;
            size_t avail = samples_in_buffer_ - min(samples_in_buffer_, N_samples_to_drop_);
            if (avail)
            {
                out_samples_gen_accum_ += avail;
                audio_chunk * out = insert_chunk(avail*channel_count_);
                out->set_data(out_buffer_, avail, channel_count_, out_rate_, channel_map_);
                samples_in_buffer_ -= avail; /* samples_in_buffer_ = N_samples_to_drop_ */
                as_memmove2(out_buffer_, 0, out_buffer_, avail, samples_in_buffer_);
            }
        }
        close();
        return;
    }
}

double dsp_rate::get_latency()
{
    if (sample_rate_ && out_rate_)
        return double(in_samples_accum_)/double(sample_rate_) - double(out_samples_gen_accum_)/double(out_rate_);
    else return 0;
    // warning: sample_rate_ and out_rate_ can be from previous track (?),
    // but in_samples_accum_ == out_samples_gen_accum_ == 0, so it's ok.
}
