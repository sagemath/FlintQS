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

// =====================================================================
// INTERFACE:
//
// void ChineseInit(void)
//         - Initialise variables for chinese
//
// void modmul(mpz_t ab, mpz_t a, mpz_t b, mpz_t p)
//         - sets ab to a*b modulo p
//
// void chinese(mpz_t res, mpz_t n, mpz_t x1, mpz_t x2, mpz_t n1, mpz_t n2)
//         - sets n to n1*n2
//         - sets res mod n to a value congruent to x1 mod n1 and x2 mod n2
//
// ======================================================================
 
extern void ChineseInit(void);
extern void modmul(mpz_t, mpz_t, mpz_t, mpz_t);
extern void chinese(mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);
 


