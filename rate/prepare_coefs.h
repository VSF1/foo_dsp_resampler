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

#define coef(coef_p, interp_order, fir_len, phase_num, coef_interp_num, fir_coef_num) coef_p[(fir_len) * ((interp_order) + 1) * (phase_num) + ((interp_order) + 1) * (fir_coef_num) + (interp_order - coef_interp_num)]

static sox_sample_t * prepare_coefs(raw_coef_t const * coefs, int num_coefs,
    int num_phases, int interp_order, int multiplier)
{
  int i, j, length = num_coefs4 * num_phases;
  sox_sample_t * result = lsx_aligned_calloc(length * (interp_order + 1), sizeof(*result));
  double fm1 = coefs[0], f1 = 0, f2 = 0;

  for (i = num_coefs - 1; i >= 0; --i)
    for (j = num_phases - 1; j >= 0; --j) {
      double f0 = fm1, b = 0, c = 0, d = 0; /* = 0 to kill compiler warning */
      int pos = i * num_phases + j - 1;
      fm1 = pos > 0 ? coefs[pos - 1] * multiplier : 0;
      switch (interp_order) {
        case 1: b = f1 - f0; break;
        case 2: b = f1 - (.5 * (f2+f0) - f1) - f0; c = .5 * (f2+f0) - f1; break;
        case 3: c=.5*(f1+fm1)-f0;d=(1/6.)*(f2-f1+fm1-f0-4*c);b=f1-f0-d-c; break;
        default: if (interp_order) assert(0);
      }
      #define coef_coef(x) \
        coef(result, interp_order, num_coefs4, j, x, num_coefs4 - 1 - i)
      coef_coef(0) = f0;
      if (interp_order > 0) coef_coef(1) = b;
      if (interp_order > 1) coef_coef(2) = c;
      if (interp_order > 2) coef_coef(3) = d;
      #undef coef_coef
      f2 = f1, f1 = f0;
    }

#define c result
/* TODO: cleaner code, maybe out-of-place permute */
/*  t[j*n + k] = a[off + j + k*n]; then a[off + j + k*n] = t[j + k*n]; */
#if defined(SSE3_)
  if (interp_order==1) {
    sox_sample_t t;
    for(i=0; i < num_coefs4*num_phases/2; i++) {
      t = c[4*i + 1]; c[4*i + 1] = c[4*i + 2]; c[4*i + 2] = t;
    }
  }
  else if (interp_order==2) {
    sox_sample_t t;
    for(i=0; i < num_coefs4*num_phases/2; i++) {
      t = c[6*i + 1]; c[6*i + 1] = c[6*i + 3]; c[6*i + 3] = c[6*i + 4]; c[6*i + 4] = c[6*i + 2]; c[6*i + 2] = t;
    }
  }
  else if (interp_order==3) {
    sox_sample_t t;
    for(i=0; i < num_coefs4*num_phases/2; i++) {
      t = c[8*i + 1]; c[8*i + 1] = c[8*i + 4]; c[8*i + 4] = c[8*i + 2]; c[8*i + 2] = t;
      t = c[8*i + 3]; c[8*i + 3] = c[8*i + 5]; c[8*i + 5] = c[8*i + 6]; c[8*i + 6] = t;
    }
  }
#elif defined(SSE_)
  if (interp_order==1) {
    sox_sample_t t;
    for(i=0; i < num_coefs4*num_phases/4; i++) {
      t = c[8*i + 1]; c[8*i + 1] = c[8*i + 2]; c[8*i + 2] = c[8*i + 4]; c[8*i + 4] = t;
      t = c[8*i + 3]; c[8*i + 3] = c[8*i + 6]; c[8*i + 6] = c[8*i + 5]; c[8*i + 5] = t;
    }
  }
  else if (interp_order==2) {
    sox_sample_t t;
    for(i=0; i < num_coefs4*num_phases/4; i++) {
      t = c[12*i + 1]; c[12*i + 1] = c[12*i + 3]; c[12*i + 3] = c[12*i + 9]; c[12*i + 9] = c[12*i + 5 ]; c[12*i + 5 ] = c[12*i + 4]; c[12*i + 4] = t;
      t = c[12*i + 2]; c[12*i + 2] = c[12*i + 6]; c[12*i + 6] = c[12*i + 7]; c[12*i + 7] = c[12*i + 10]; c[12*i + 10] = c[12*i + 8]; c[12*i + 8] = t;
    }
  }
  else if (interp_order==3) {
    sox_sample_t t;
    for(i=0; i < num_coefs4*num_phases/4; i++) {
      t = c[16*i + 1 ]; c[16*i + 1 ] = c[16*i + 4 ]; c[16*i + 4 ] = t;
      t = c[16*i + 2 ]; c[16*i + 2 ] = c[16*i + 8 ]; c[16*i + 8 ] = t;
      t = c[16*i + 3 ]; c[16*i + 3 ] = c[16*i + 12]; c[16*i + 12] = t;
      t = c[16*i + 6 ]; c[16*i + 6 ] = c[16*i + 9 ]; c[16*i + 9 ] = t;
      t = c[16*i + 7 ]; c[16*i + 7 ] = c[16*i + 13]; c[16*i + 13] = t;
      t = c[16*i + 11]; c[16*i + 11] = c[16*i + 14]; c[16*i + 14] = t;
    }
  }
#endif
#undef c

  return result;
}

#undef coef
