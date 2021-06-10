/* libSoX Library Public Interface
 *
 * Copyright 1999-2012 Chris Bagwell and SoX Contributors.
 *
 * This source code is freely redistributable and may be used for
 * any purpose.  This copyright notice must be maintained.
 * Chris Bagwell And SoX Contributors are not responsible for
 * the consequences of using this software.
 */

#ifndef SOX_H
#define SOX_H

#include <stddef.h> /* Ensure NULL etc. are available throughout SoX */
#include <stdlib.h>

#include <stdint.h>

typedef double sox_rate_t;

/* Boolean type, assignment (but not necessarily binary) compatible with C++ bool */
typedef enum sox_bool {
    sox_bool_dummy = -1, /* Ensure a signed type */
    sox_false = 0, /**< False = 0. */
    sox_true = 1   /**< True = 1. */
} sox_bool;

#endif
