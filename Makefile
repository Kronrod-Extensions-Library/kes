GMP_LIB_DIR=/userdata/raoulb/lib/lib
GMP_INCLUDE_DIR=/userdata/raoulb/lib/include
MPFR_LIB_DIR=/userdata/raoulb/lib/lib
MPFR_INCLUDE_DIR=/userdata/raoulb/lib/include
FLINT_LIB_DIR=/userdata/raoulb/lib/lib
FLINT_INCLUDE_DIR=/userdata/raoulb/lib/include/flint
#FLINT_CPIMPORT=\"/userdata/raoulb/lib/share/flint/CPimport.txt\"
ARB_LIB_DIR=/userdata/raoulb/lib/lib
ARB_INCLUDE_DIR=/userdata/raoulb/lib/include


POLY ?= HERMITE
DIMENSION ?= 1
PRINTLOG ?= 1
CFG= -D${POLY} -DDIMENSION=${DIMENSION} -DPRINTLOG=${PRINTLOG}


CC=gcc
CFLAGS=-std=c11 -Wall -Werror -pedantic -O2 -funroll-loops -fopenmp -mpopcnt


CPP=g++
CPPFLAGS=-std=c++11 -Wall -Werror -pedantic -O2 -funroll-loops -fopenmp -mpopcnt -fpermissive -Wno-sign-compare


INC=-I$(CURDIR) -I$(GMP_INCLUDE_DIR) -I$(MPFR_INCLUDE_DIR) -I$(FLINT_INCLUDE_DIR) -I$(ARB_INCLUDE_DIR)

LIB=-L$(CURDIR) -L$(ARB_LIB_DIR) -L$(FLINT_LIB_DIR) -L$(GMP_LIB_DIR) -L$(MPFR_LIB_DIR) -larb -lflint -lgmp -lmpfr -lpthread -lm


all: kes ekes rekes quadrature genzkeister test enumtest

quadrature: quadrature.c *.h
	$(CC) $(CFLAGS) $(CFG) $(INC) libkes2.h quadrature.c $(LIB) -o quadrature

kes: kes.c *.h
	$(CC) $(CFLAGS) $(CFG) $(INC) libkes2.h kes.c $(LIB) -o kes

ekes: kes_enumerate.c *.h
	$(CC) $(CFLAGS) $(CFG) $(INC) libkes2.h kes_enumerate.c $(LIB) -o ekes

rekes: kes_rec_enumerate.c *.h
	$(CC) $(CFLAGS) $(CFG) $(INC) libkes2.h kes_rec_enumerate.c $(LIB) -o rekes

genzkeister: genzkeister.cpp *.h
	$(CPP) $(CPPFLAGS) $(CFG) $(INC) genzkeister.h genzkeister.cpp $(LIB) -o genzkeister

test: test.c *.h
	$(CC) $(CFLAGS) $(INC) libkes2.h test.c $(LIB) -o test

enumtest: enumtest.cpp enumerators.h
	$(CPP) $(CPPFLAGS) $(CFG) $(INC) enumerators.h enumtest.cpp $(LIB) -o enumtest

clean:
	rm kes ekes rekes quadrature genzkeister test enumtest
