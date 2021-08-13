/*============================================================================
    Copyright 2006 William Hart    

    This file is part of SIMPQS.

    SIMPQS is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    SIMPQS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMPQS; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

============================================================================*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "F2matrix.h"

static u_int32_t bitPattern[]  =
{
  0x80000000, 0x40000000, 0x20000000, 0x10000000,
  0x08000000, 0x04000000, 0x02000000, 0x01000000,
  0x00800000, 0x00400000, 0x00200000, 0x00100000,
  0x00080000, 0x00040000, 0x00020000, 0x00010000,
  0x00008000, 0x00004000, 0x00002000, 0x00001000,
  0x00000800, 0x00000400, 0x00000200, 0x00000100,
  0x00000080, 0x00000040, 0x00000020, 0x00000010,
  0x00000008, 0x00000004, 0x00000002, 0x00000001
};

void insertEntry(matrix m, u_int32_t i, u_int32_t j)
{
     m[i][j / 32] |= bitPattern[j % 32];
     
     return;
}

void xorEntry(matrix m, u_int32_t i, u_int32_t j)
{
     m[i][j / 32] ^= bitPattern[j % 32];
     
     return;
}

u_int32_t getEntry(matrix m, u_int32_t i, u_int32_t j)
{
     return m[i][j / 32] & bitPattern[j % 32];
}

void swapRows(matrix m, u_int32_t x, u_int32_t y)
{
     row temp;
     
     temp = m[x];
     m[x] = m[y];
     m[y] = temp;
     
     return;
}


void clearRow(matrix m, u_int32_t numcols, u_int32_t row)
{
    int32_t dwords = numcols/32;
    
    if (numcols%32) dwords++;
    memset(m[row],0,dwords*4);
    
    return; 
}

void displayRow(matrix m, u_int32_t row, u_int32_t numPrimes)
{
     int32_t length = numPrimes/32;
     if (numPrimes%32) length++;
     length*=64;
     
     printf("[");
     for (int32_t j = 0; j < length/2; j++)
     {
        if (getEntry(m,row,j)) printf("1");
        else printf("0");
     }
     printf("  ");
     for (int32_t j = length/2; j < length; j++)
     {
        if (getEntry(m,row,j)) printf("1");
        else printf("0");
     }
     printf("]\n");
     return;
}

void xorRows(matrix m, u_int32_t source, u_int32_t dest, u_int32_t length)
{
  u_int32_t i, q, r;
  row x = m[dest];
  row y = m[source];
  
  r = length % 8; q = length - r;
  for (i=0; i < q; i += 8)
  {
    x[i] ^= y[i]; x[1+i] ^= y[1+i]; x[2+i] ^= y[2+i]; x[3+i] ^= y[3+i];
    x[4+i] ^= y[4+i]; x[5+i] ^= y[5+i]; x[6+i] ^= y[6+i]; x[7+i] ^= y[7+i];
  }
  switch (r)
  {
    case 7: x[i] ^= y[i]; i++;
    case 6: x[i] ^= y[i]; i++;
    case 5: x[i] ^= y[i]; i++;
    case 4: x[i] ^= y[i]; i++;
    case 3: x[i] ^= y[i]; i++;
    case 2: x[i] ^= y[i]; i++;
    case 1: x[i] ^= y[i]; i++;
  }
  
  return;
}

matrix constructMat(u_int32_t cols, u_int32_t rows)
{
     static matrix m;
     
     u_int32_t dwords = cols/32;
     if (cols%32) dwords++;
     m = (row *) calloc(sizeof(row),rows);
     if (m==NULL) 
     {
       printf("Unable to allocate memory for matrix!\n");
       exit(1);
     }

     for (u_int32_t i = 0; i < rows; i++)
     {
         m[i] = (row) calloc(2*dwords,sizeof(u_int32_t));  //two matrices, side by side
     }
     if (m[rows-1]==NULL) 
     {
        printf("Unable to allocate memory for matrix!\n");
        exit(1);
     }
     
     for (u_int32_t i = 0; i < rows; i++)  //make second matrix identity, i.e. 1's along diagonal
     {
         insertEntry(m,i,i+32*dwords);
     }
     
     return m;
}

//===========================================================================
// gaussReduce:
//
// Function: Apply Gaussian elimination to a matrix.
//
//===========================================================================
u_int32_t gaussReduce(matrix m, u_int32_t numPrimes, u_int32_t relSought,int32_t extras)
{
     static u_int32_t rowUpto = 0;
     static u_int32_t irow;
     static u_int32_t length = (numPrimes+extras)/32;
     
     if (numPrimes%32) length++;
     length*=2;
     
     for (int32_t icol = numPrimes-1; icol >= 0; icol--)
     {
         irow = rowUpto;
         while ((irow < relSought)&&(getEntry(m,irow,icol)==0UL)) irow++;
         if (irow < relSought) 
         {
             
             swapRows(m,rowUpto,irow);
             
             for (u_int32_t checkRow = rowUpto+1; checkRow<relSought; checkRow++)
             {
                 if (getEntry(m,checkRow,icol)!=0UL) 
                 {
                    xorRows(m,rowUpto,checkRow,length);
                 }
             }
             
             rowUpto++;
         }
     }
          
     return rowUpto;
}

