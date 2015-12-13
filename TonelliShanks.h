/*============================================================================
    Copyright 2006 William Hart    

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

============================================================================*/

#include <stdlib.h>

// =====================================================================
// INTERFACE:
//
// void TonelliInit(void) 
//         - Initialises variables
//
// int sqrtmod(mpz_t asqrt, mpz_t a, mpz_t p) 
//         - Tonelli-Shanks: sets asqrt to a square root of a modulo p
//         - Return: 0 if a is not a square mod p, 1 otherwise.
//
// void sqrtmodpk(mpz_t res, mpz_t z, mpz_t a, mpz_t p, int k)
//         - Given a square root z, of a mod p (from Tonelli-Shanks say)
//         - sets res to a square root of a mod p^k
//
// ========================================================================

extern void TonelliInit(void);
extern int32_t sqrtmod(mpz_t, mpz_t, mpz_t);
extern inline void sqrtmodpow(mpz_t, mpz_t, mpz_t, mpz_t);
extern void sqrtmodpk(mpz_t, mpz_t, mpz_t, mpz_t, int32_t);

 
