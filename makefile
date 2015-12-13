#============================================================================
#    Copyright 2006 William Hart    
#
#    This file is part of FLINT.
#
#    FLINT is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    FLINT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with FLINT; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
#
#============================================================================

CPP  = g++
OBJ  = lprels.o ModuloArith.o TonelliShanks.o F2matrix.o lanczos.o QS.o $(RES)
LINKOBJ  = lprels.o ModuloArith.o TonelliShanks.o F2matrix.o lanczos.o QS.o $(RES)
#LIBS =  -L"home/dmharvey/gmp/install/lib" -lgmp 
#CXXINCS = -I"home/dmharvey/gmp/install/include" 
LIBS = -L"/usr/lib" -lgmp
CXXINCS = -I"/usr/include"
BIN  = QuadraticSieve QuadraticSieve.exe
CXXFLAGS = $(CXXINCS) -Wall -Wno-sign-compare -fomit-frame-pointer -O2
CXXFLAGS2 = $(CXXINCS) -Wall -Wno-sign-compare -fomit-frame-pointer -O3
#CXXFLAGS = $(CXXINCS) -ansi -Wall -Wno-sign-compare -march=athlon-xp -fomit-frame-pointer -O2
#CXXFLAGS2 = $(CXXINCS) -ansi -Wall -Wno-sign-compare -march=athlon-xp -fomit-frame-pointer -O3
#CXXFLAGS = $(CXXINCS) -Wall -Wno-sign-compare -march=opteron -fomit-frame-pointer -O2 
#CXXFLAGS2 = $(CXXINCS) -Wall -Wno-sign-compare -march=opteron -fomit-frame-pointer -O3
RM = rm -f

.PHONY: all clean clean-custom

all: QuadraticSieve

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) -ansi $(LINKOBJ) -o "QuadraticSieve" $(LIBS)

ModuloArith.o: ModuloArith.cpp
	$(CPP) -ansi -c ModuloArith.cpp -o ModuloArith.o $(CXXFLAGS)

TonelliShanks.o: TonelliShanks.cpp
	$(CPP) -ansi -c TonelliShanks.cpp -o TonelliShanks.o $(CXXFLAGS)

F2matrix.o: F2matrix.cpp
	$(CPP) -ansi -c F2matrix.cpp -o F2matrix.o $(CXXFLAGS2)
	
lanczos.o: lanczos.c
	$(CPP) -ansi -c lanczos.c -o lanczos.o $(CXXFLAGS2)       

lprels.o: lprels.c
	$(CPP) -ansi -c lprels.c -o lprels.o $(CXXFLAGS2)

QS.o: QS.cpp
	$(CPP) -ansi -c QS.cpp -o QS.o $(CXXFLAGS)
