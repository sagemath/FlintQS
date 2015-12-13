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

//  -------------------------------------------------------
//
//  ModuloArith.cpp
//
//  Provides Functions for doing Modulo Arithmetic
//
//  -------------------------------------------------------
 
#include <gmp.h>
#include "ModuloArith.h"

mpz_t restemp;  //chinese variables
mpz_t ntemp;     
mpz_t chtemp;  

void modmul(mpz_t ab, mpz_t a, mpz_t b, mpz_t p)
{
     mpz_mul(ab,a,b);
     mpz_fdiv_r(ab,ab,p);
}

void ChineseInit(void)
{
    mpz_init(restemp);
    mpz_init(ntemp);
    mpz_init(chtemp);
    
    return;
}
    
void chinese(mpz_t res, mpz_t n, mpz_t x1, mpz_t x2, mpz_t n1, mpz_t n2)
{
     mpz_mul(ntemp,n1,n2);
     mpz_invert(restemp,n2,n1);
     modmul(restemp,restemp,n2,ntemp);
     modmul(restemp,restemp,x1,ntemp);
     
     mpz_invert(chtemp,n1,n2);
     modmul(chtemp,chtemp,n1,ntemp);
     modmul(chtemp,chtemp,x2,ntemp);
     
     mpz_add(res,restemp,chtemp);
     mpz_fdiv_r(res,res,ntemp);
     mpz_set(n,ntemp);
     
     return;
}

