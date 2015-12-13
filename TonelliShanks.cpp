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

// -------------------------------------------------------
//
// TonelliShanks.cpp
//
// Provides Tonelli-Shanks square root mod p, and mod p^k
//
// -------------------------------------------------------
 
#include <gmp.h>
#include "TonelliShanks.h"
#include "ModuloArith.h"  //need multiplication mod p

mpz_t two; //variables for sqrtmod
mpz_t p1;
mpz_t b;
mpz_t g;
mpz_t xsq;
mpz_t mk;
mpz_t bpow;
mpz_t gpow;
     
mpz_t inv;  //variables for sqrtmodpow
mpz_t tempsqpow;
              
mpz_t pk;  //variable for sqrtmodpk

void TonelliInit(void)
{
    mpz_init(two);
    mpz_init(p1);
    mpz_init(b);
    mpz_init(g);
    mpz_init(xsq);
    mpz_init(mk);
    mpz_init(bpow);
    mpz_init(gpow);
     
    mpz_init(inv);
    mpz_init(tempsqpow);
              
    mpz_init(pk);
    
    return;
}

int32_t sqrtmod(mpz_t asqrt, mpz_t a, mpz_t p) 
{
     int32_t r,k,m;
     
     if (mpz_kronecker(a,p)!=1) 
     {
         mpz_set_ui(asqrt,0);
         return 0;   //return 0 if a is not a square mod p
     }
     
     mpz_set_ui(two,2);
     
     mpz_sub_ui(p1,p,1);
     r = mpz_remove(p1,p1,two);
     mpz_powm(b,a,p1,p);
     for (k=2; ;k++)
     {
         if (mpz_ui_kronecker(k,p) == -1) break;
     }
     mpz_set_ui(mk,k);
     mpz_powm(g,mk,p1,p);
     mpz_add_ui(p1,p1,1);
     mpz_divexact_ui(p1,p1,2);
     mpz_powm(xsq,a,p1,p);
     if (!mpz_cmp_ui(b,1)) 
     {
          mpz_set(asqrt,xsq);
          
          return 1;
     }
     
     while (mpz_cmp_ui(b,1))
     {
           mpz_set(bpow,b);
           for (m=1; (m<=r-1) && (mpz_cmp_ui(bpow,1));m++)
           {
               mpz_powm_ui(bpow,bpow,2,p);
           }
           mpz_set(gpow,g);
           for (int32_t i = 1;i<= r-m-1;i++)
           {
               mpz_powm_ui(gpow,gpow,2,p);
           };
           modmul(xsq,xsq,gpow,p);
           mpz_powm_ui(gpow,gpow,2,p);
           modmul(b,b,gpow,p);
           mpz_set(gpow,g);
           r = m;
     }
          
     mpz_set(asqrt,xsq);
     
     return 1;
}

inline void sqrtmodpow(mpz_t res, mpz_t z, mpz_t a, mpz_t pk)
{
     mpz_mul_ui(inv,z,2);
     mpz_invert(inv,inv,pk);
     mpz_set(tempsqpow,a);
     mpz_submul(tempsqpow,z,z);
     mpz_fdiv_r(tempsqpow,tempsqpow,pk);
     modmul(tempsqpow,tempsqpow,inv,pk);
     mpz_add(tempsqpow,tempsqpow,z);
     mpz_set(res,tempsqpow);
     
     return;
} 

void sqrtmodpk(mpz_t res, mpz_t z, mpz_t a, mpz_t p, int32_t k)
{
     mpz_set(res,z);
     mpz_set(pk,p);
     for (int32_t i = 2;i<=k;i++)
     {
            mpz_mul(pk,pk,p);
            sqrtmodpow(res,res,a,pk);
     }
     
     return;
}
