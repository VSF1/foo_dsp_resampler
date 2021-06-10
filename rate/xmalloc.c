/* SoX Memory allocation functions
 *
 * Copyright (c) 2005-2006 Reuben Thomas.  All rights reserved.
 * Copyright (C) 2008-10 lvqcl
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

#include "xmalloc.h"
#undef malloc
#undef calloc
#undef realloc
#undef free

#include <stdlib.h>
#include "rate_i.h"

void (*malloc_error_handler)(void) = NULL;

void set_malloc_error_handler(void (*h)(void))
{
    malloc_error_handler = h;
}


void * lsx_malloc(size_t size)
{
    void * p = malloc(size);
    if (p == NULL) malloc_error_handler();
    return p;
}

void *lsx_realloc(void *ptr, size_t newsize)
{
    void * p;
    if (newsize == 0) { lsx_free(ptr); return NULL; }
    if (ptr == NULL) return lsx_malloc(newsize);

    if ((p = realloc(ptr, newsize)) == NULL) {
        /* lsx_free(ptr);  -- we don't want double free */
        malloc_error_handler();
    }
    return p;
}

void* lsx_aligned_malloc(size_t size)
{
    void * p = _aligned_malloc(size, LSX_ALIGN);
    if (p == NULL) malloc_error_handler();
    return p;
}

void * lsx_aligned_realloc(void * ptr, size_t newsize, size_t oldsize)
{
    void * p;
    if (newsize == 0) { lsx_aligned_free(ptr); return NULL; }
    if (ptr == NULL) return lsx_aligned_malloc(newsize);

    p = _aligned_realloc(ptr, newsize, LSX_ALIGN);
    if (p) return p;

    if ((p = _aligned_malloc(newsize, LSX_ALIGN)) == NULL) {
        /* lsx_aligned_free(ptr);  -- we don't want double free */
        malloc_error_handler();
    }

    memcpy(p, ptr, oldsize);
    lsx_aligned_free(ptr);
    return p;
}
