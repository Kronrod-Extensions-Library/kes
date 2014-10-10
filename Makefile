GMP_LIB_DIR=/userdata/raoulb/lib/lib
GMP_INCLUDE_DIR=/userdata/raoulb/lib/include
MPFR_LIB_DIR=/userdata/raoulb/lib/lib
MPFR_INCLUDE_DIR=/userdata/raoulb/lib/include
FLINT_LIB_DIR=/userdata/raoulb/lib/lib
FLINT_INCLUDE_DIR=/userdata/raoulb/lib/include/flint
ARB_LIB_DIR=/userdata/raoulb/lib/lib
ARB_INCLUDE_DIR=/userdata/raoulb/lib/include


POLY ?= HERMITE
PRINTLOG ?= 1
CFG=-DPRINTLOG=${PRINTLOG} -D${POLY}


CC=gcc
CFLAGS=-std=c11 -fopenmp -pedantic -Wall -Werror -O2 -funroll-loops -mpopcnt -DFLINT_CPIMPORT=\"/userdata/raoulb/lib/share/flint/CPimport.txt\"


CPP=g++
CPPFLAGS=-std=c++11 -Wall -O2 -funroll-loops -mpopcnt -pedantic -fpermissive -Wno-sign-compare -Wno-unused-variable -DFLINT_CPIMPORT=\"/userdata/raoulb/lib/share/flint/CPimport.txt\"


INCS=-I$(CURDIR) -I$(GMP_INCLUDE_DIR) -I$(MPFR_INCLUDE_DIR) -I$(FLINT_INCLUDE_DIR) -I$(ARB_INCLUDE_DIR)

LIBS=-L$(CURDIR) -L$(ARB_LIB_DIR) -L$(FLINT_LIB_DIR) -L$(GMP_LIB_DIR) -L$(MPFR_LIB_DIR) -larb -lflint -lmpfr -lgmp -lpthread -lm


bkes:
	$(CC) $(CFLAGS) $(CFG) $(INCS) libkes2.h kes.c $(LIBS) -o kes

ekes:
	$(CC) $(CFLAGS) $(CFG) $(INCS) libkes2.h kes_enumerate.c $(LIBS) -o kes_enumerate

rekes:
	$(CC) $(CFLAGS) $(CFG) $(INCS) libkes2.h kes_rec_enumerate.c $(LIBS) -o kes_rec_enumerate

quad:
	$(CC) $(CFLAGS) $(CFG) $(INCS) libkes2.h quadrature.c $(LIBS) -o quadrature

gk:
	$(CPP) $(CPPFLAGS) $(CFG) $(INCS) genzkeister.cpp $(LIBS) -o gkq

test:
	$(CC) $(CFLAGS) $(INCS) libkes2.h test.c $(LIBS) -o test

clean:
	rm kes kes_enumerate kes_rec_enumerate quadrature test
