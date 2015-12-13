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

Description:

This is a relatively fast implementation of the self-initialising quadratic sieve.
If you manage to improve the code, the author would like to hear about it.

Contact: hart_wb {at-thingy} yahoo.com
================================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include <sys/times.h>
#include <limits.h>

#include "TonelliShanks.h"
#include "ModuloArith.h"
#include "F2matrix.h"
#include "lanczos.h"
#include "lprels.h"

//===========================================================================
//Uncomment these for various pieces of debugging information

#define COUNT    // Shows the number of relations generated and curves used during sieving
//#define RELPRINT     // Shows the actual factorizations of the relations
//#define ERRORS   // Error if relation should be divisible by a prime but isn't 
//#define POLS     // Shows the polynomials being used by the sieve
//#define ADETAILS // Prints some details about the factors of the A coefficients of the polys
//#define LARGESTP // Prints the size of the largest factorbase prime
//#define CURPARTS // Prints the number of curves used and number of partial relations
//#define TIMING //displays some relative timings, if feature is available
//#define REPORT //report sieve size, multiplier and number of primes used

//===========================================================================
//Architecture dependent fudge factors

#if ULONG_MAX == 4294967295U
#define SIEVEMASK 0xC0C0C0C0U
#define MIDPRIME 1500
#define SIEVEDIV 1
#elif ULONG_MAX == 18446744073709551615U
#define SIEVEMASK 0xC0C0C0C0C0C0C0C0U
#define MIDPRIME       1500 
#define SIEVEDIV 1 
#endif

#define CACHEBLOCKSIZE 64000 //Should be a little less than the L1/L2 cache size
                             //and a multiple of 64000
#define MEDIUMPRIME    900   
#define SECONDPRIME    6000 //This should be lower for slower machines
#define FUDGE          0.15 //Every program needs a mysterious fudge factor

#define MINDIG 40 //Will not factor numbers with less than this number of decimal digits

#define PREFETCH(addr,n) __builtin_prefetch((unsigned long*)addr+n,0,1)

//===========================================================================
//Knuth-Schroeppel multipliers and a macro to count them

static const unsigned long multipliers[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 
                                                23, 29, 31, 37, 41, 43};

#define NUMMULTS (sizeof(multipliers)/sizeof(unsigned long))

//===========================================================================
// Large prime cutoffs
                   
unsigned long largeprimes[] = 
{
     250000, 300000, 370000, 440000, 510000, 580000, 650000, 720000, 790000, 8600000, //40-49
     930000, 1000000, 1700000, 2400000, 3100000, 3800000, 4500000, 5200000, 5900000, 6600000, //50-59
     7300000, 8000000, 8900000, 10000000, 11300000, 12800000, 14500000, 16300000, 18100000, 20000000, //60-69
     22000000, 24000000, 27000000, 32000000, 39000000,  //70-74
     53000000, 65000000, 75000000, 87000000, 100000000, //75-79
     114000000, 130000000, 150000000, 172000000, 195000000, //80-84
     220000000, 250000000, 300000000, 350000000, 400000000, //85-89
     450000000, 500000000 //90-91
};

//============================================================================
// Number of primes to use in factor base, given the number of decimal digits specified
unsigned long primesNo[] = 
{
     1500, 1500, 1600, 1700, 1750, 1800, 1900, 2000, 2050, 2100, //40-49
     2150, 2200, 2250, 2300, 2400, 2500, 2600, 2700, 2800, 2900, //50-59
     3000, 3150, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, //60-69 
     9500, 10000, 11500, 13000, 15000, //70-74
     17000, 24000, 27000, 30000, 37000, //75-79
     45000, 47000, 53000, 57000, 58000,  //80-84
     59000, 60000, 64000, 68000, 72000,  //85-89
     76000, 80000 //90-91
};

//============================================================================
// First prime actually sieved for
unsigned long firstPrimes[] = 
{
     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, //40-49
     9, 8, 9, 9, 9, 9, 10, 10, 10, 10, //50-59
     10, 10, 11, 11, 12, 12, 13, 14, 15, 17, //60-69  //10
     19, 21, 22, 22, 23, //70-74
     24, 25, 25, 26, 26, //75-79
     27, 27, 27, 27, 28, //80-84
     28, 28, 28, 29, 29, //85-89
     29, 29 //90-91
};

//============================================================================
// Logs of primes are rounded and errors accumulate; this specifies how great an error to allow
unsigned long errorAmounts[] = 
{
     16, 17, 17, 18, 18, 19, 19, 19, 20, 20, //40-49
     21, 21, 21, 22, 22, 22, 23, 23, 23, 24, //50-59
     24, 24, 25, 25, 25, 25, 26, 26, 26, 26, //60-69 //24
     27, 27, 28, 28, 29, //70-74
     29, 30, 30, 30, 31, //75-79
     31, 31, 31, 32, 32, //80-84
     32, 32, 32, 33, 33, //85-89
     33, 33 //90-91
};

//============================================================================
// This is the threshold the sieve value must exceed in order to be considered for smoothness
unsigned long thresholds[] = 
{
     66, 67, 67, 68, 68, 68, 69, 69, 69, 69, //40-49
     70, 70, 70, 71, 71, 71, 72, 72, 73, 73, //50-59
     74, 74, 75, 75, 76, 76, 77, 77, 78, 79, //60-69 //74
     80, 81, 82, 83, 84, //70-74
     85, 86, 87, 88, 89, //75-79
     91, 92, 93, 93, 94, //80-84 
     95, 96, 97, 98, 100, //85-89
     101, 102 //90-91  
};

//============================================================================
// Size of sieve to use divided by 2, given the number of decimal digits specified
//N.B: probably optimal if chosen to be a multiple of 32000, though other sizes are supported
unsigned long sieveSize[] = 
{
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //40-49
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //50-59
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //60-69
     32000, 32000, 64000, 64000, 64000, //70-74
     96000, 96000, 96000, 128000, 128000, //75-79
     160000, 160000, 160000, 160000, 160000, //80-84 
     192000, 192000, 192000, 192000, 192000, //85-89
     192000, 192000 //90-91
};

// Athlon tuning parameters
/*unsigned long sieveSize[] = 
{
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //40-49
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //50-59
     64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, //60-69
     64000, 64000, 64000, 64000, 64000, //70-74
     128000, 128000, 128000, 128000, 128000, //75-79
     160000, 160000, 160000, 160000, 160000, //80-84 
     192000, 192000, 192000, 192000, 192000, //85-89
     192000, 192000 //90-91
};*/

//============================================================================
long decdigits; //number of decimal digits of n
unsigned long secondprime; //min(numprimes, SECONDPRIME) = cutoff for using flags when sieving
unsigned long firstprime;  //first prime actually sieved with
unsigned char errorbits;  //first prime actually sieved with
unsigned char threshold;  //sieve threshold cutoff for smooth relations
unsigned long midprime;
unsigned long largeprime;

unsigned long * factorBase; //array of factor base primes
unsigned long numPrimes; //number of primes in factor base
unsigned long relSought; //number of relations sought, i.e. a "few" more than numPrimes
unsigned char * primeSizes; //array of sizes in bits, of the factor base primes
unsigned char * sieve; //actual array where sieving takes place
unsigned char * * offsets; //offsets for each prime to use in sieve 
unsigned char * * offsets2; //offsets for each prime to use in sieve (we switch between these)
unsigned long relsFound =0; //number of relations found so far
unsigned long potrels = 0; //potential relations (including duplicates)
unsigned char * flags; //flags used for speeding up sieving for large primes
unsigned long partials = 0; //number of partial relations
unsigned long Mdiv2; //size of sieving interval divide 2 
unsigned long mat2off; //offset of second square block in matrix

mpz_t * sqrts; //square roots of n modulo each prime in the factor base

mpz_t n; //number to be factored 
mpz_t res; //smooth values which are trial factored

mpz_t temp, temp2, temp3; //temporary variables
mpz_t q,r; //quotient and remainder

//Variables used for keeping time

unsigned long clockstart;
unsigned long clocktotal = 0;  

//Variables used by the modular inversion macro function
long u1, u3;
long v1, v3;
long t1, t3, quot;

//Variable used for random function
unsigned long randval = 2994439072U;

//==========================================================================================
//Timing: provides some relative timings on X86 machines running gcc

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))

#ifdef TIMING
#define TIMES
#endif

double counterfirst[4];
double countertotal[4] = {0,0,0,0};

static unsigned counthi = 0;
static unsigned countlo = 0;

void counterasm(unsigned *hi, unsigned *lo)
{
 asm("rdtsc; movl %%edx,%0; movl %%eax,%1" 
 : "=r" (*hi), "=r" (*lo) 
 : 
 : "%edx", "%eax");
}

double getcounter()
{
   double total;

   counterasm(&counthi, &countlo);

   total = (double) counthi * (1 << 30) * 4 + countlo;
   return total;
}

#endif
   
/*========================================================================
   Modular Inversion:

   Function: GMP has a modular inverse function, but believe it or not, 
             this clumsy implementation is apparently quite a bit faster. 
             It inverts the value a, modulo the prime p, using the extended 
             gcd algorithm.

========================================================================*/

inline unsigned long modinverse(unsigned long a, unsigned long p)
{
   u1=1; u3=a;
   v1=0; v3=p;
   t1=0; t3=0;
   while (v3)
   {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            { 
               t1 = u1; u1 = v1; v1 = t1;
               t3 = u3; u3 = v3; v3 = t3;
            } else 
            {
               t1 = u1 - v1; u1 = v1; v1 = t1;
               t3 = u3 - v3; u3 = v3; v3 = t3;
            }
         } else if (quot < (v3<<1))
         {  
            t1 = u1 - (v1<<1); u1 = v1; v1 = t1;
            t3 = u3 - (v3<<1); u3 = v3; v3 = t3;
         } else
         {
            t1 = u1 - v1*3; u1 = v1; v1 = t1;
            t3 = u3 - v3*3; u3 = v3; v3 = t3;
         }
      } else
      {
         quot=u3/v3;
         t1 = u1 - v1*quot; u1 = v1; v1 = t1;
         t3 = u3 - v3*quot; u3 = v3; v3 = t3;
      }
   } 
   
   if (u1<0) u1+=p;
   
   return u1;
}

/*=========================================================================
   Knuth_Schroeppel Multiplier:
 
   Function: Find the best multiplier to use (allows 2 as a multiplier).
             The general idea is to find a multiplier k such that kn will
             be faster to factor. This is achieved by making kn a square 
             modulo lots of small primes. These primes will then be factor
             base primes, and the more small factor base primes, the faster
             relations will accumulate, since they hit the sieving interval
             more often. 
 
==========================================================================*/
unsigned long knuthSchroeppel(mpz_t n)
{
    float bestFactor = -10.0f;
    unsigned long multiplier = 1;
    unsigned long nmod8;
    float factors[NUMMULTS];
    float logpdivp;
    mpz_t prime, r, mult;
    long kron, multindex;
    
    mpz_init(prime);
    mpz_init(r);
    mpz_init(mult);
    
    nmod8 = mpz_fdiv_r_ui(r,n,8);
    
    for (multindex = 0; multindex < NUMMULTS; multindex++)
    {
       long mod = nmod8*multipliers[multindex]%8;
       factors[multindex] = 0.34657359; // ln2/2 
       if (mod == 1) factors[multindex] *= 4.0;   
       if (mod == 5) factors[multindex] *= 2.0;   
       factors[multindex] -= (log((float) multipliers[multindex]) / 2.0);
    }
    
    mpz_set_ui(prime,3);
    while (mpz_cmp_ui(prime,10000)<0)
    {
          logpdivp = log((float)mpz_get_ui(prime)) / mpz_get_ui(prime);
          kron = mpz_kronecker(n,prime);
          for (multindex = 0; multindex < NUMMULTS; multindex++)
          {
              mpz_set_ui(mult,multipliers[multindex]);
              switch (kron*mpz_kronecker(mult,prime))
              {
                 case 0:
                 {
                      factors[multindex] += logpdivp;
                 } break;
                 case 1:
                 {
                      factors[multindex] += 2.0*logpdivp;
                 } break;
                 default: break;
              }
          }
          
          mpz_nextprime(prime,prime);
    }
    
    for (multindex=0; multindex<NUMMULTS; multindex++)
    {
      if (factors[multindex] > bestFactor)
      { 
        bestFactor = factors[multindex];
        multiplier = multipliers[multindex];
      }
    } 
    
    mpz_clear(prime);
    mpz_clear(r);
    mpz_clear(mult);
    
    return multiplier;
}



/*========================================================================
   Initialize Quadratic Sieve:
  
   Function: Initialises the global gmp variables.

========================================================================*/
void initSieve(void)
{
    mpz_init(n);
    mpz_init(temp); 
    mpz_init(temp2);
    mpz_init(temp3);
    mpz_init(res);
    mpz_init(q);
    mpz_init(r);
    
    return;
}

/*========================================================================
   Compute Factor Base:
 
   Function: Computes primes p up to B for which n is a square mod p,  
   allocates memory and stores them in an array pointed to by factorBase
   Returns: number of primes actually in the factor base

========================================================================*/
void computeFactorBase(mpz_t n, unsigned long B,unsigned long multiplier)
{
     mpz_t currentPrime;
     unsigned long primesinbase = 0;
     
     factorBase = (unsigned long *) calloc(sizeof(unsigned long),B); 
     
     factorBase[primesinbase] = multiplier;
     primesinbase++;
     if (multiplier!=2)
     {
        factorBase[primesinbase] = 2;
        primesinbase++;
     }
     mpz_init_set_ui(currentPrime,3);
     while (primesinbase < B)
     {
          if (mpz_kronecker(n,currentPrime)==1)
          {
              factorBase[primesinbase] = mpz_get_ui(currentPrime);
              primesinbase++;
          } 
          mpz_nextprime(currentPrime,currentPrime);
     }
#ifdef LARGESTP
     gmp_printf("Largest prime less than %Zd\n",currentPrime);
#endif
      
     mpz_clear(currentPrime);
     return;
}

/*===========================================================================
   Compute Prime Sizes:
 
   Function: Computes the size in bits of each prime in the factor base
     allocates memory for an array, primeSizes, to store the sizes
     stores the size for each of the numPrimes primes in the array 
 
===========================================================================*/
void computeSizes(unsigned long numPrimes)
{
     primeSizes = (unsigned char *) calloc(sizeof(unsigned char),numPrimes);
     for (unsigned long i = 0; i<numPrimes; i++)
     {
         primeSizes[i]=(unsigned char)floor(log((double)factorBase[i])/log(2.0)-FUDGE+0.5);
     }
     
     return;
}

/*===========================================================================
   Tonelli-Shanks:

   Function: Performs Tonelli-Shanks on n mod every prime in the factor base
      allocates memory for the results to be stored in the array sqrts

===========================================================================*/
void tonelliShanks(unsigned long numPrimes,mpz_t n)
{
     sqrts = (mpz_t *) calloc(sizeof(mpz_t),numPrimes); 
     mpz_array_init(sqrts[0],numPrimes,8*sizeof(unsigned long));
     
     for (unsigned long i = 1; i<numPrimes; i++) 
     {
         mpz_set_ui(temp,factorBase[i]);
         sqrtmod(sqrts[i],n,temp);
     }
     
     return;
}

/*==========================================================================
   evaluateSieve:

   Function: searches sieve for relations and sticks them into a matrix, then
             sticks their X and Y values into two arrays XArr and YArr

===========================================================================*/
void evaluateSieve(unsigned long ** relations, unsigned long ctimesreps, unsigned long M, unsigned char * sieve, mpz_t A, mpz_t B, mpz_t C, unsigned long * soln1, unsigned long * soln2, long polyadd, unsigned long * polycorr, mpz_t * XArr, unsigned long * aind, long min, long s,unsigned long multiplier, long * exponents, la_col_t* colarray, unsigned long * factors, char * rel_str, FILE* LPNEW,FILE* RELS)
{
     long i,j;
     register unsigned long k;
     unsigned long exponent, vv;
     unsigned char extra;
     register unsigned long modp;
     unsigned long * sieve2;
     unsigned char bits;
     long numfactors;
     unsigned long factnum;
     char * last_ptr;
     char Q_str[200];
     char X_str[200];
     
     i = 0;
     j=0;
     sieve2 = (unsigned long *) sieve;
#ifdef POLS
     gmp_printf("%Zdx^2%+Zdx\n%+Zd\n",A,B,C);
#endif
     
     while (j<M/sizeof(unsigned long))
     {
        do
        {
           while (!(sieve2[j] & SIEVEMASK)) j++;
           i=j*sizeof(unsigned long);
           j++;
           while ((i<j*sizeof(unsigned long))&&(sieve[i] < threshold)) i++;
        } while (sieve[i] < threshold);
           
        if (i<M) 
        {
           mpz_set_ui(temp,i+ctimesreps);
           mpz_sub_ui(temp,temp,Mdiv2); //X
              
           mpz_set(temp3,B);  //B
           mpz_addmul(temp3,A,temp);  //AX+B
           mpz_add(temp2,temp3,B);  //AX+2B
           mpz_mul(temp2,temp2,temp);  //AX^2+2BX
           mpz_add(res,temp2,C);  //AX^2+2BX+C
              
           bits=mpz_sizeinbase(res,2);
           bits-=errorbits;
              
           numfactors=0;
              
           extra = 0;
           if (factorBase[0]!=1)
           {
              mpz_set_ui(temp,factorBase[0]);
              exponent = mpz_remove(res,res,temp);
              exponents[0] = exponent;
              if (exponent) 
              { 
                 extra+=primeSizes[0];
              }
           }
             
           mpz_set_ui(temp,factorBase[1]);
           exponent = mpz_remove(res,res,temp);
           exponents[1] = exponent;
           extra+=exponent;
                
           for (k = 2; k<firstprime; k++)
           {
              modp=(i+ctimesreps)%factorBase[k];
                
              if (soln2[k]!=0xFFFFFFFFl)
              {
                 if ((modp==soln1[k]) || (modp==soln2[k]))
                 {
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
             
#ifdef ERRORS
                    if (exponent==0) printf("Error!\n");
#endif
                    extra+=primeSizes[k];
#ifdef RELPRINT
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%ld",exponent);
#endif
                    exponents[k] = exponent;
                 } else exponents[k] = 0;
              } else
              {
                 mpz_set_ui(temp,factorBase[k]);
                 exponent = mpz_remove(res,res,temp);
                 if (exponent) extra+=primeSizes[k];
#ifdef RELPRINT
                 if (exponent > 0) gmp_printf(" %Zd",factorBase[k]);
                 if (exponent > 1) printf("^%ld",exponent);
#endif
                 exponents[k] = exponent;
              }  
           }  
           factnum = 0;
           sieve[i]+=extra;
           if (sieve[i] >= bits)
           {
              vv=((unsigned char)1<<(i&7));
              for (k = firstprime; (k<secondprime)&&(extra<sieve[i]); k++)
              {
                 modp=(i+ctimesreps)%factorBase[k];
                 if (soln2[k]!=0xFFFFFFFFl)
                 {
                    if ((modp==soln1[k]) || (modp==soln2[k]))
                    {
                       extra+=primeSizes[k];
                       mpz_set_ui(temp,factorBase[k]);
                       exponent = mpz_remove(res,res,temp);
              
#ifdef ERRORS
                       if (exponent==0) printf("Error!\n");
#endif
                        
#ifdef RELPRINT
                       if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                       if (exponent > 1) printf("^%ld",exponent);
#endif
                       factors[factnum+1] = k;
                       factors[factnum] = exponent;
                       factnum+=2;
                    }  
                 } else
                 {
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
                    if (exponent) extra+=primeSizes[k];
                        
#ifdef RELPRINT
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%ld",exponent);
#endif
                    if (exponent) 
                    { 
                       factors[factnum+1] = k;
                       factors[factnum] = exponent;
                       factnum+=2;
                    }
                 }  
              }  
              
              for (k = secondprime; (k<numPrimes)&&(extra<sieve[i]); k++)
              {
                 if (flags[k]&vv)
                 {
                    modp=(i+ctimesreps)%factorBase[k];
                    if ((modp==soln1[k]) || (modp==soln2[k]))
                    {
                       mpz_set_ui(temp,factorBase[k]);
                       exponent = mpz_remove(res,res,temp);
#ifdef ERRORS
                       if (exponent==0) printf("Error!\n");
#endif
                       extra+=primeSizes[k];
#ifdef RELPRINT
                       if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                       if (exponent > 1) printf("^%ld",exponent);
#endif
                       factors[factnum+1] = k;
                       factors[factnum] = exponent;
                       factnum+=2;
                    }  
                 }  
              }  
              
              last_ptr = rel_str;
              if (mpz_cmp_ui(res,1000)>0)
              {
                 if (mpz_cmp_ui(res,largeprime)<0) 
                 {
                    for (unsigned long i = 0; i < firstprime; i++)
                    {
                        if (exponents[i]) add_factor(&last_ptr, (unsigned long) exponents[i], (unsigned long) i);
                    }
                    for (unsigned long i = 0; i < factnum; i+=2)
                    {
                        add_factor(&last_ptr, (unsigned long) factors[i], (unsigned long) factors[i+1]);
                    }
                    for (long i =0; i<s; i++)
                    {
                        add_factor(&last_ptr, (unsigned long) 1, (unsigned long) aind[i]+min);
                    }
                    
                    add_0(&last_ptr);
                    gmp_sprintf(X_str, "%Zd\0", temp3);
                    gmp_sprintf(Q_str, "%Zd\0", res);
                    fprintf(LPNEW, "%s @ %s :%s\n", Q_str, X_str, rel_str);
                    partials++;
                 }
#ifdef RELPRINT
                 gmp_printf(" %Zd\n",res);
#endif
              } else
              { 
                 mpz_neg(res,res);
                 if (mpz_cmp_ui(res,1000)>0)
                 {
                    if (mpz_cmp_ui(res,largeprime)<0) 
                    {
                       for (unsigned long i = 0; i < firstprime; i++)
                       {
                          if (exponents[i]) add_factor(&last_ptr, (unsigned long) exponents[i], (unsigned long) i);
                       }
                       for (unsigned long i = 0; i < factnum; i+=2)
                       {
                          add_factor(&last_ptr, (unsigned long) factors[i], (unsigned long) factors[i+1]);
                       }
                       for (long i =0; i<s; i++)
                       {
                          add_factor(&last_ptr, (unsigned long) 1, (unsigned long) aind[i]+min);
                       }
                    
                       add_0(&last_ptr);
                       gmp_sprintf(X_str, "%Zd\0", temp3);
                       gmp_sprintf(Q_str, "%Zd\0", res);
                       fprintf(LPNEW, "%s @ %s :%s\n", Q_str, X_str, rel_str);
                    
                       partials++;
                    }
#ifdef RELPRINT
                    gmp_printf(" %Zd\n",res);
#endif
                 } else 
                 {
#ifdef RELPRINT
                    printf("....R\n");
#endif
                    for (long i = 0; i<firstprime; i++) 
                    {
                       if (exponents[i]) add_factor(&last_ptr, (unsigned long) exponents[i], (unsigned long) i);
                    }
                    for (unsigned long i = 0; i < factnum; i+=2)
                    {
                       add_factor(&last_ptr, (unsigned long) factors[i], (unsigned long) factors[i+1]);
                    }
                    for (long i =0; i<s; i++)
                    {
                       add_factor(&last_ptr, (unsigned long) 1, (unsigned long) aind[i]+min);
                    }
                       
                    add_0(&last_ptr);
                    gmp_sprintf(X_str, "%Zd\0", temp3);
                    fprintf(RELS, "%s :%s\n", X_str, rel_str);
                                           
                    potrels++;
                 } 
              }  
           } else 
           {
#ifdef RELPRINT
              printf("\r                                                                    \r");
#endif
                 
           }   
           i++;
              
        };      
     }  
     
     return;
}


/*=============================================================================
   Sieve:

   Function: Allocates space for a sieve of M integers and sieves the interval
             starting at start

=============================================================================*/
void sieveInterval(unsigned long M, unsigned long numPrimes, unsigned char * sieve, long last, long first, long polyadd, unsigned long * soln1, unsigned long * soln2, unsigned long * polycorr, unsigned char * * offsets, unsigned char * * offsets2)
{
     register unsigned char currentprimesize; 
     register unsigned long currentprime;
     unsigned char * position2;
     register unsigned char * position;
     register long diff;
     unsigned char * end;
     unsigned long ptimes4;
     long correction;
     
     end = sieve+M;
     
     if (first)
     {
        for (unsigned long prime=1; prime<firstprime; prime++) 
        {
            if (soln2[prime] == 0xFFFFFFFF) continue;
            currentprime = factorBase[prime];
            correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
            soln1[prime]+=correction;
            while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
            soln2[prime]+=correction;
            while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
        }  
     }
     
     for (unsigned long prime=firstprime; prime<MEDIUMPRIME; prime++) 
     {
        if (soln2[prime] == 0xFFFFFFFF) continue;
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        if (first)
        {
           correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
           soln1[prime]+=correction;
           while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
           soln2[prime]+=correction;
           while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;

           position = sieve+soln1[prime];
           position2 = sieve+soln2[prime];
        } else
        {
           position = offsets[prime];
           position2 = offsets2[prime];
        }
        diff=position2-position;
        
        ptimes4 = currentprime*4;
        register unsigned char * bound=end-ptimes4;
        while (bound - position > 0)  
        {  
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
        }
        while ((end - position > 0)&&(end - position - diff > 0))
        { 
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
              
        }
        position2 = position+diff;
        if (end - position2 > 0)
        { 
              (* position2)+=currentprimesize, position2+=currentprime;
        }
        if (end - position > 0)
        { 
              (* position)+=currentprimesize, position+=currentprime;
        }
        
        if (!last)
        {
           offsets[prime] = position;
           offsets2[prime] = position2;
        }
     } 
     
      for (unsigned long prime=MEDIUMPRIME; prime<midprime; prime++) 
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        if (first)
        {
           correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
           soln1[prime]+=correction;
           while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
           soln2[prime]+=correction;
           while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
           
           position = sieve+soln1[prime];
           position2 = sieve+soln2[prime];
        } else
        {
           position = offsets[prime];
           position2 = offsets2[prime];
        }
        diff=position2-position;
           
        ptimes4 = 2*currentprime;
        register unsigned char * bound=end-ptimes4;
        while (bound - position > 0)  
        {  
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
        }
        position2 = position+diff;
        while ((end - position > 0)&&(end - position2 > 0))
        { 
              (* position)+=currentprimesize, position+=currentprime, (* position2)+=currentprimesize, position2+=currentprime;
              
        }
        
        if (end - position2 > 0)
        { 
              (* position2)+=currentprimesize, position2+=currentprime;
        }
        if (end - position > 0)
        { 
              (* position)+=currentprimesize, position+=currentprime;
        }
        if (!last)
        {
           offsets[prime] = position;
           offsets2[prime] = position2;
        }
     } 
     
     return;
}

/*===========================================================================
   Sieve 2:

   Function: Second sieve for larger primes

=========================================================================== */
void sieve2(unsigned long M, unsigned long numPrimes, unsigned char * sieve, long last, long first, long polyadd, unsigned long * soln1, unsigned long * soln2, unsigned long * polycorr, unsigned char * * offsets, unsigned char * * offsets2)
{
     register unsigned char currentprimesize; 
     register unsigned long currentprime;
     register unsigned char * position2;
     register unsigned char * position;
     unsigned char * end;
     long correction;
     
     memset(sieve,0,M*sizeof(unsigned char));
     memset(flags,0,numPrimes*sizeof(unsigned char));
     end = sieve+M;
     *end = 255; //sentinel to speed up sieve evaluators inner loop

     for (unsigned long prime=midprime; prime<secondprime; prime++) 
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
        soln1[prime]+=correction;
        while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
        soln2[prime]+=correction;
        while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
           
        position = sieve+soln1[prime];
        position2 = sieve+soln2[prime];
           
        while ((end - position > 0)&&(end - position2 > 0))
        { 
             (* position)+=currentprimesize, position+=currentprime, (* position2)+=currentprimesize, position2+=currentprime;
        }
        
        if (end - position2 > 0)
        { 
              (* position2)+=currentprimesize;
        }
        if (end - position > 0)
        { 
              (* position)+=currentprimesize;
        }        
     }
     
     for (unsigned long prime=secondprime; prime<numPrimes; prime++) 
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
        soln1[prime]+=correction;
        while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
        soln2[prime]+=correction;
        while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
           
        position = sieve+soln1[prime];
        position2 = sieve+soln2[prime];
           
        while (end - position > 0)
        { 
              flags[prime]|=((unsigned char)1<<((position-sieve)&7)), (* position)+=currentprimesize, position+=currentprime;
        }
        
        while (end - position2 > 0)
        { 
              flags[prime]|=((unsigned char)1<<((position2-sieve)&7)), (* position2)+=currentprimesize, position2+=currentprime;
        }
     }
     
     return;
}

/*============================================================================

   random: 
           
   Function: Generates a pseudo-random integer between 0 and n-1 inclusive
   
============================================================================*/
unsigned long random(unsigned long upto)
{
   randval = ((u_int64_t)randval*1025416097U+286824428U)%(u_int64_t)4294967291U;
   return randval%upto;
}


/*============================================================================
   mainRoutine:

   Function: Generates the polynomials, initialises and calls the sieve, 
             implementing cache blocking (breaking the sieve interval into
             small blocks for the small primes.

============================================================================*/
void mainRoutine(unsigned long Mdiv2, mpz_t n, unsigned long multiplier)
{
    mpz_t A; mpz_init(A);
    mpz_t B; mpz_init(B);
    mpz_t C; mpz_init(C);
    mpz_t D; mpz_init(D);
    mpz_t temp; mpz_init(temp);
    mpz_t temp2; mpz_init(temp2);
    mpz_t q; mpz_init(q);
    mpz_t r; mpz_init(r);
    mpz_t Bdivp2; mpz_init(Bdivp2);
    mpz_t factor; mpz_init(factor);
          
    unsigned long u1;
     
    long s, fact, span, min;
    unsigned long p;     
    unsigned long reps;
     
    unsigned long curves = 0; 
     
    unsigned long ** relations;
    long * primecount;

    long * exponents = (long *) calloc(firstprime,sizeof(long));
    if (exponents==NULL) 
    {
       printf("Unable to allocate memory!\n");
       abort();
    }
    unsigned long factors[200];
    char rel_str[MPQS_STRING_LENGTH];
    
    unsigned long totcomb = 0;
    
    unsigned long next_cutoff = (relSought - 1)/40 +1;
    unsigned long next_inc = next_cutoff;
    
    FILE * LPNEW;
    FILE * LPRELS;
    FILE * COMB;
    FILE * FNEW;
    FILE * RELS;
    FILE * FRELS;
    FILE * FLPRELS;
    LPNEW = flint_fopen("lpnew","w");
    LPRELS = flint_fopen("lprels","w");
    RELS = flint_fopen("rels","w");
    FNEW = flint_fopen("fnew","w");
    fclose(FNEW);
    FLPRELS = flint_fopen("flprels","w");
    fclose(FLPRELS);
    FRELS = flint_fopen("frels","w");
    fclose(LPRELS);
    fclose(FRELS);
    
 
#ifdef TIMES
    counterfirst[2] = getcounter();
#endif
    s = mpz_sizeinbase(n,2)/28+1;
     
    unsigned long * aind = (unsigned long*) calloc(sizeof(unsigned long),s);  
    unsigned long * amodp = (unsigned long*) calloc(sizeof(unsigned long),s); 
    unsigned long * Ainv = (unsigned long*) calloc(sizeof(unsigned long),numPrimes); 
    unsigned long * soln1 = (unsigned long*) calloc(sizeof(unsigned long),numPrimes); 
    unsigned long * soln2 = (unsigned long*) calloc(sizeof(unsigned long),numPrimes); 
    unsigned long ** Ainv2B = (unsigned long**) calloc(sizeof(unsigned long*),s);
    if (Ainv2B==NULL) 
    {
       printf("Unable to allocate memory!\n");
       abort();
    }
    for (long i=0; i<s; i++)
    {
       Ainv2B[i] = (unsigned long *) calloc(sizeof(unsigned long),numPrimes);
       if (Ainv2B[i]==NULL) 
       {
          printf("Unable to allocate memory!\n");
          abort();
       }
    } 
     
    mpz_t * Bterms = (mpz_t *)calloc(sizeof(mpz_t),s);
    mpz_array_init(*Bterms,s,mpz_sizeinbase(n,2)); 
    
    la_col_t* colarray = (la_col_t*)calloc(relSought,sizeof(la_col_t)); //initialise a zeroed array of column structures
    
    relsFound = 0;
     
    mpz_t XArr[relSought];
    for(unsigned long i = 0; i < relSought; i++) mpz_init(XArr[i]);
     
    sieve = (unsigned char *) calloc(sizeof(unsigned char),Mdiv2*2+4); //one dword extra for sentinel to speed up sieve evaluation loop
    if (sieve==NULL) 
    {
       printf("Unable to allocate memory for sieve!\n");
       abort();
    }                
 
     
    flags = (unsigned char*) calloc(sizeof(unsigned char),numPrimes);
     
    offsets = (unsigned char * *)calloc(numPrimes,sizeof(unsigned char *));
    offsets2 = (unsigned char * *)calloc(numPrimes,sizeof(unsigned char *));
     
    relations = (unsigned long * *) calloc(relSought,sizeof(unsigned long *));
    for (long i = 0; i < relSought; i++)
    {
       relations[i] = (unsigned long *) calloc(200, sizeof(unsigned long));
    }
     
    primecount = (long *) calloc(numPrimes, sizeof(long));
 
//Compute min A_prime and A_span
     
    mpz_mul_ui(temp,n,2);
    mpz_sqrt(temp,temp);
    mpz_div_ui(temp2,temp,Mdiv2);
    mpz_root(temp,temp2,s);
    for (fact = 0; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++); 
    span = numPrimes/s/s/2;
    min=fact-span/2;
    while ((fact*fact)/min - min < span) {min--;}
    
#ifdef ADETAILS
    printf("s = %ld, fact = %ld, min = %ld, span = %ld\n",s,fact,min,span);
#endif
     
//Compute first polynomial and adjustments
    
    while (relsFound + totcomb < relSought)
    {
        long i,j;
        long ran;
           
        mpz_set_ui(A,1);
        
        for (i = 0; i < s-1; )
        {
           j=-1L;
           ran = span/2+random(span/2);
           while (j!=i)
           {
              ran++;
              for (j=0;((j<i)&&(aind[j]!=ran));j++);
           }
           aind[i] = ran;
           mpz_mul_ui(A,A,factorBase[ran+min]);
           i++;
           if (i < s-1)
           {
              j=-1L;
              ran = ((min+span/2)*(min+span/2))/(ran+min) - random(10)-min;
              while (j!=i)
              {
                 ran++;
                 for (j=0;((j<i)&&(aind[j]!=ran));j++);
              }
              aind[i] = ran;
              mpz_mul_ui(A,A,factorBase[ran+min]);
              i++;
           } 
        } 
        
        mpz_div(temp,temp2,A);
        for (fact = 1; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++); 
        fact-=min;
        do
        {
           for (j=0;((j<i)&&(aind[j]!=fact));j++);
           fact++;
        } while (j!=i);
        fact--;
        aind[i] = fact;
        mpz_mul_ui(A,A,factorBase[fact+min]); 
           
        for (long i=0; i<s; i++)
        {
           p = factorBase[aind[i]+min];
           mpz_div_ui(temp,A,p);
	       amodp[i] = mpz_fdiv_r_ui(temp,temp,p);
                          
	       mpz_set_ui(temp,modinverse(mpz_get_ui(temp),p));
           mpz_mul(temp,temp,sqrts[aind[i]+min]);
           mpz_fdiv_r_ui(temp,temp,p);
           if (mpz_cmp_ui(temp,p/2)>0)
           {
              mpz_sub_ui(temp,temp,p);
              mpz_neg(temp,temp);
           }
           mpz_mul(temp,temp,A);
           mpz_div_ui(Bterms[i],temp,p);     
        }
            
        mpz_set(B,Bterms[0]);
        for (long i=1; i<s; i++)
        {
           mpz_add(B,B,Bterms[i]);
        }
            
        for (long i=0; i<numPrimes; i++)
        {
           p = factorBase[i];
	       Ainv[i] = modinverse(mpz_fdiv_r_ui(temp,A,p),p);
             
           for (long j=0; j<s; j++)
           {
              mpz_fdiv_r_ui(temp,Bterms[j],p);
              mpz_mul_ui(temp,temp,2*Ainv[i]);
              Ainv2B[j][i] = mpz_fdiv_r_ui(temp,temp,p);
           }
             
           mpz_fdiv_r_ui(temp,B,p);
           mpz_sub(temp,sqrts[i],temp);
           mpz_add_ui(temp,temp,p);
           mpz_mul_ui(temp,temp,Ainv[i]);
           mpz_add_ui(temp,temp,Mdiv2);
           soln1[i] = mpz_fdiv_r_ui(temp,temp,p);
           mpz_sub_ui(temp,sqrts[i],p);
           mpz_neg(temp,temp);
           mpz_mul_ui(temp,temp,2*Ainv[i]);
           soln2[i] = mpz_fdiv_r_ui(temp,temp,p)+soln1[i];
        }  
            
        for (long polyindex=1; polyindex<(1<<(s-1))-1; polyindex++)
        {
           long j,polyadd;
           unsigned long * polycorr;
           for (j=0; j<s; j++)
           {
              if (((polyindex>>j)&1)!=0) break;
           }
           if ((polyadd = (((polyindex>>j)&2)!=0)))
           {
              mpz_add(B,B,Bterms[j]);
              mpz_add(B,B,Bterms[j]);
           } else
           {
              mpz_sub(B,B,Bterms[j]); 
              mpz_sub(B,B,Bterms[j]); 
           }
           polycorr = Ainv2B[j];
             
           long index;
           for (long j=0; j<s; j++)
           {
              index = aind[j]+min;
              p = factorBase[index];
              mpz_fdiv_r_ui(D,n,p*p);
              mpz_fdiv_r_ui(Bdivp2,B,p*p);
              mpz_mul_ui(temp,Bdivp2,amodp[j]);
              mpz_realloc2(temp3,64);
	          mpz_fdiv_r_ui(temp,temp,p);
	          u1 = modinverse(mpz_fdiv_r_ui(temp,temp,p),p);        
              mpz_mul(temp,Bdivp2,Bdivp2);
              mpz_sub(temp,temp,D);
              mpz_neg(temp,temp);
              mpz_div_ui(temp,temp,p);
              mpz_mul_ui(temp,temp,u1);
              mpz_add_ui(temp,temp,Mdiv2);
              mpz_add_ui(temp,temp,p);
	          soln1[index]=mpz_fdiv_r_ui(temp,temp,p);
              soln2[index]=0xFFFFFFFFl;
           }
           
// Count the number of polynomial curves used so far and compute the C coefficient of our polynomial
     
           curves++;
             
           mpz_mul(C,B,B);
           mpz_sub(C,C,n);
           mpz_divexact(C,C,A);
           
// Do the sieving and relation collection
     
           mpz_set_ui(temp,Mdiv2*2);
           mpz_fdiv_qr_ui(q,r,temp,CACHEBLOCKSIZE);
#ifdef TIMES
           counterfirst[3] = getcounter();
#endif
           sieve2(mpz_get_ui(temp),numPrimes,sieve,1,1,polyadd,soln1,soln2,polycorr,offsets,offsets2);
#ifdef TIMES
           countertotal[3]+=(getcounter()-counterfirst[3]);
           counterfirst[0] = getcounter();
#endif
           sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve,0,1,polyadd,soln1,soln2,polycorr,offsets,offsets2);
           if (mpz_cmp_ui(q,1)>0)
           {
              for (reps = 1;reps < mpz_get_ui(q)-1; reps++)
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,0,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              }
              if (mpz_cmp_ui(r,0)==0)
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,1,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              } else
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,0,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
                 reps++;
                 sieveInterval(mpz_get_ui(r),numPrimes,sieve+CACHEBLOCKSIZE*reps,1,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              }
           }
           
#ifdef TIMES  
           countertotal[0]+=(getcounter()-counterfirst[0]);
           counterfirst[1] = getcounter();
#endif
           evaluateSieve(relations,0,mpz_get_ui(temp),sieve,A,B,C,soln1, soln2, polyadd, polycorr,XArr,aind,min,s,multiplier,exponents,colarray,factors,rel_str,LPNEW,RELS);
           
           if (2*potrels >= next_cutoff)
           {
              fclose(LPNEW); 
              sort_lp_file("lpnew");
              COMB = flint_fopen("comb","w");
              mergesort_lp_file("lprels", "lpnew", "tmp", COMB);
              fclose(COMB);
              LPNEW = flint_fopen("lpnew","w");
              
              fclose(RELS);
              sort_lp_file("rels");
              relsFound = mergesort_lp_file("frels","rels","tmp2",NULL);
              RELS = flint_fopen("rels","w");
              
              COMB = flint_fopen("comb", "r");
              FNEW = flint_fopen("fnew","w");
              combine_large_primes(numPrimes, COMB, FNEW, n, factor);
              fclose(FNEW);
              fclose(COMB);
              sort_lp_file("fnew");
              totcomb = mergesort_lp_file("flprels","fnew","tmp3",NULL);
#ifdef COUNT
              printf("%ld full relations, %ld combined relations\n",relsFound,totcomb);
#endif
              if ((next_cutoff < relSought) && (next_cutoff + next_inc/2 >= relSought)) 
                next_inc = next_inc/2;
              next_cutoff += next_inc;
           }
#ifdef TIMES
           countertotal[1]+=(getcounter()-counterfirst[1]);
#endif
       }  
     
#ifdef COUNT
        if (curves%20==0) printf("%ld curves.\n",(long)curves);
#endif
    }
     
#ifdef CURPARTS
    printf("%ld curves, %ld partials.\n",(long)curves,(long)partials);
#endif
             
#ifdef REPORT
    printf("Done with sieving!\n");
#endif
     
   unsigned long ncols = relSought;
   unsigned long nrows = numPrimes;

#ifdef ERRORS
   for (unsigned long j = relsFound; j<relSought; j++)
   {
       if (colarray[j].weight != 0) printf("Dirty at col %ld\n",j);
   }
#endif

   fclose(LPNEW); 
   fclose(RELS);
   
#ifdef COUNT
   printf("%ld relations found in total!\n",totcomb+relsFound);
#endif

   FRELS = flint_fopen("frels","r");
   relsFound = 0;
   read_matrix(relations, FRELS, colarray, &relsFound, relSought, XArr, n, factorBase);
   fclose(FRELS);
   FLPRELS = flint_fopen("flprels","r");
   read_matrix(relations, FLPRELS, colarray, &relsFound, relSought, XArr, n, factorBase);
   fclose(FLPRELS);
   
#ifdef ERRORS   
   for (unsigned long j = 0; j<relSought; j++)
   {
      if (colarray[j].orig != j) 
      {
         printf("Column numbering error, %ld\n",j);
         colarray[j].orig = j;
      }
   }
   
   for (unsigned long j = 0; j<relSought; j++)
      for (unsigned long i = 0; i<colarray[j].weight; i++) 
        if (colarray[j].data[i] > numPrimes) printf("Error prime too large: %ld\n",colarray[j].data[i]);

   mpz_t test1;
   mpz_init(test1);
   mpz_t test2;
   mpz_init(test2);
   mpz_t test3;
   mpz_init(test3);
   unsigned long * exps = (unsigned long *) malloc(numPrimes*sizeof(unsigned long));
   for (unsigned long j = 0; j<relSought; j++)
   {
      for (unsigned long i = 0; i < numPrimes; i++) exps[i] = 0;
      mpz_set_ui(test1,1);
      for (unsigned long i = 1; i<=relations[j][0]; i++)
      {
        mpz_mul_ui(test1,test1,factorBase[relations[j][i]]);
        exps[relations[j][i]]++;
      }
      mpz_mod(test1,test1,n);
      mpz_mul(test2,XArr[j],XArr[j]);
      mpz_mod(test2,test2,n);
      if (mpz_cmp(test1,test2)!=0)
      {
         mpz_add(test3,test1,test2);
         if (mpz_cmp(test3,n)!=0) 
         {
            gmp_printf("%Zd !=\n%Zd\nin column %ld and\n",test3,n,j);
            gmp_printf("%Zd !=\n%Zd\n\n",test1,test2);
         }
      }
      for (unsigned long i = 0; i < colarray[j].weight; i++) 
      {
         if (exps[colarray[j].data[i]]%2 != 1) printf("Col %ld, row %ld incorrect\n",j,i);
         exps[colarray[j].data[i]] = 0;
      }
      for (unsigned long i = 0; i < numPrimes; i++) if (exps[i]%2==1) printf("exps[%ld] is not even in row %ld\n",i,j);
   }
   free(exps);
   mpz_clear(test1);
   mpz_clear(test2);
   mpz_clear(test3);
#endif
 
   reduce_matrix(&nrows, &ncols, colarray);
 
#ifdef ERRORS  
   exps = (unsigned long *) malloc(numPrimes*sizeof(unsigned long));
   for (unsigned long j = 0; j<ncols; j++)
   {
       for (unsigned long i = 0; i < numPrimes; i++) exps[i] = 0;
       for (unsigned long i = 1; i<=relations[colarray[j].orig][0]; i++)
       {
          exps[relations[colarray[j].orig][i]]++;
       }
       for (unsigned long i = 0; i < colarray[j].weight; i++) 
       {
          for (unsigned long k = 0; k < i; k++)
             if (colarray[j].data[i] == colarray[j].data[k]) printf("Duplicate in column %ld: %ld, %ld\n",j,i,k);
          if (exps[colarray[j].data[i]]%2 != 1) printf("Col %ld, row %ld incorrect\n",j,i);
          exps[colarray[j].data[i]] = 0;
       }
       for (unsigned long i = 0; i < numPrimes; i++) if (exps[i]%2==1) printf("exps[%ld] is not even in row %ld\n",i,j);
       
   }
#endif
     
     u_int64_t* nullrows;
     do {
         nullrows = block_lanczos(nrows, 0, ncols, colarray);
     } while (nullrows == NULL);
     
     long i,j, mask;
     
     for (i = 0, mask = 0; i < ncols; i++)
		mask |= nullrows[i];

	for (i = j = 0; i < 64; i++) {
		if (mask & ((u_int64_t)(1) << i))
			j++;
	}
#ifdef REPORT
	printf ("%ld nullspace vectors found.\n",j);
#endif

 
#ifdef ERRORS  
   exps = (unsigned long *) malloc(numPrimes*sizeof(unsigned long));
   for (unsigned long j = 0; j<ncols; j++)
   {
       for (unsigned long i = 0; i < numPrimes; i++) exps[i] = 0;
       for (unsigned long i = 1; i<=relations[colarray[j].orig][0]; i++)
       {
          exps[relations[colarray[j].orig][i]]++;
       }
       for (unsigned long i = 0; i < colarray[j].weight; i++) 
       {
          if ((exps[colarray[j].data[i]]%2) != 1) printf("Col %ld, row %ld incorrect\n",j,i);
          exps[colarray[j].data[i]] = 0;
       }
       for (unsigned long i = 0; i < numPrimes; i++) if ((exps[i]%2)==1) printf("exps[%ld] is not even in row %ld\n",i,j);
   }
#endif
    

#ifdef TIMES
    countertotal[2]+=(getcounter()-counterfirst[2]);
    printf("Total time = %.0f\n",countertotal[2]);
    printf("Polynomial generation time = %.0f\n",countertotal[2]-countertotal[0]-countertotal[1]-countertotal[3]);
    printf("Small prime sieve time = %.0f\n",countertotal[0]);
    printf("Large prime sieve time = %.0f\n",countertotal[3]);
    printf("Evaluate sieve time = %.0f\n",countertotal[1]); 
#endif
    
// We want factors of n, not kn, so divide out by the multiplier
     
    mpz_div_ui(n,n,multiplier);
    
// Now do the "square root" and GCD steps hopefully obtaining factors of n
    printf("FACTORS:\n");
    for (long l = 0;l<64;l++)
    {
        while (!(mask & ((u_int64_t)(1) << l))) l++;
        mpz_set_ui(temp,1);
        mpz_set_ui(temp2,1);
        memset(primecount,0,numPrimes*sizeof(unsigned long));
        for (long i = 0; i<ncols; i++)
        {
           if (getNullEntry(nullrows,i,l)) 
           {
              mpz_mul(temp2,temp2,XArr[colarray[i].orig]);
              for (long j=1; j<=relations[colarray[i].orig][0]; j++)
              {
                 primecount[relations[colarray[i].orig][j]]++;
              }
           }
           if (i%30==0) mpz_mod(temp2,temp2,n);
        }
        for (long j = 0; j<numPrimes; j++) 
        {
           mpz_set_ui(temp3,factorBase[j]);
           mpz_pow_ui(temp3,temp3,primecount[j]/2);
           mpz_mul(temp,temp,temp3);
           if (j%30==0) mpz_mod(temp,temp,n);
        }
        mpz_sub(temp,temp2,temp);
        mpz_gcd(temp,temp,n);
        if ((mpz_cmp(temp,n)!=0)&&(mpz_cmp_ui(temp,1)!=0)) //print only nontrivial factors
        {
           gmp_printf("%Zd\n",temp);
        }
    }
    mpz_clear(factor);
    return;
}

/*===========================================================================
   Main Program:

   Function: Factors a user specified number using a quadratic sieve

===========================================================================*/
int main(int argc, char *argv[])
{
    unsigned long multiplier;

    initSieve(); 
    
    printf("Input number to factor [ >=40 decimal digits]: "); 
    gmp_scanf("%Zd",n);getchar();
    
    decdigits = mpz_sizeinbase(n,10);
    if (decdigits < 40) 
    {
       printf("Error in input or number has too few digits.\n");
       abort();
    }
    
    multiplier = knuthSchroeppel(n);
    mpz_mul_ui(n,n,multiplier);
    
  if (decdigits<=91) 
  {
    numPrimes=primesNo[decdigits-MINDIG];
    
    Mdiv2 = sieveSize[decdigits-MINDIG]/SIEVEDIV;
    if (Mdiv2*2 < CACHEBLOCKSIZE) Mdiv2 = CACHEBLOCKSIZE/2;
    largeprime = largeprimes[decdigits-MINDIG];
    
#ifdef REPORT
    printf("Using multiplier: %ld\n",(long)multiplier);
    printf("%ld primes in factor base.\n",(long)numPrimes);
    printf("Sieving interval M = %ld\n",(long)Mdiv2*2);
    printf("Large prime cutoff = factorBase[%ld]\n",largeprime);
#endif
    
    if (numPrimes < SECONDPRIME) secondprime = numPrimes;
    else secondprime = SECONDPRIME;
    if (numPrimes < MIDPRIME) midprime = numPrimes;
    else midprime = MIDPRIME;
    
    firstprime = firstPrimes[decdigits-MINDIG];
    errorbits = errorAmounts[decdigits-MINDIG];
    threshold = thresholds[decdigits-MINDIG];
    
  } else //all bets are off
  {
     numPrimes = 64000;
     Mdiv2 = 192000/SIEVEDIV;
     largeprime = numPrimes*10*decdigits;
     
#ifdef REPORT
    printf("Using multiplier: %ld\n",(long)multiplier);
    printf("%ld primes in factor base.\n",(long)numPrimes);
    printf("Sieving interval M = %ld\n",(long)Mdiv2*2);
    printf("Large prime cutoff = factorBase[%ld]\n",largeprime);
#endif
    
    secondprime = SECONDPRIME;
    midprime = MIDPRIME;
    firstprime = 30;
    errorbits = decdigits/4 + 2;
    threshold = 43+(7*decdigits)/10;
    
  }
    relSought = numPrimes+64; 
    computeFactorBase(n, numPrimes, multiplier);

    computeSizes(numPrimes);
    
    TonelliInit();
    tonelliShanks(numPrimes,n);
    
    mainRoutine(Mdiv2, n,multiplier);
    
    getchar();
#if defined(WINCE) || defined(macintosh)
    char * tmp_dir = NULL;
#else
    char * tmp_dir = getenv("TMPDIR");
#endif
    if (tmp_dir == NULL) tmp_dir = "./";
    char * delfile;
    
    delfile = get_filename(tmp_dir,unique_filename("comb"));
    remove(delfile);
    delfile = get_filename(tmp_dir,unique_filename("frels"));
    remove(delfile);
    delfile = get_filename(tmp_dir,unique_filename("flprels"));
    remove(delfile);
    delfile = get_filename(tmp_dir,unique_filename("lpnew"));
    remove(delfile);
    delfile = get_filename(tmp_dir,unique_filename("rels"));
    remove(delfile);
    delfile = get_filename(tmp_dir,unique_filename("fnew"));
    remove(delfile);
    delfile = get_filename(tmp_dir,unique_filename("lprels"));
    remove(delfile);
      
    return 0;
}

