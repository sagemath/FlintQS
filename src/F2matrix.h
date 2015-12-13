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
#ifndef F2MATRIX_H
#define F2MATRIX_H

#include <stdlib.h>
#include "lanczos.h"

typedef u_int32_t * row; //row of an F2 matrix
typedef row * matrix; //matrix as a list of pointers to rows

extern void insertEntry(matrix, u_int32_t, u_int32_t);
extern void xorEntry(matrix, u_int32_t, u_int32_t);
extern u_int32_t getEntry(matrix, u_int32_t, u_int32_t);
extern matrix constructMat(u_int32_t, u_int32_t);
extern void xorRows(row, row, u_int32_t);
extern void clearRow(matrix, u_int32_t, u_int32_t);
extern void swapRows(row, row);
extern u_int32_t gaussReduce(matrix, u_int32_t, u_int32_t, int32_t);
extern void displayRow(matrix, u_int32_t, u_int32_t);

#endif





 
