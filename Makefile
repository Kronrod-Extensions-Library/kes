GMP_LIB_DIR=/userdata/raoulb/lib/lib
GMP_INCLUDE_DIR=/userdata/raoulb/lib/include
MPFR_LIB_DIR=/userdata/raoulb/lib/lib
MPFR_INCLUDE_DIR=/userdata/raoulb/lib/include
FLINT_LIB_DIR=/userdata/raoulb/lib/lib
FLINT_INCLUDE_DIR=/userdata/raoulb/lib/include/flint

CC=gcc

#CFLAGS=-std=c11 -pedantic -Wall -funroll-loops -g -mpopcnt -DFLINT_CPIMPORT=\"/userdata/raoulb/lib/share/flint/CPimport.txt\"
CFLAGS=-std=c11 -pedantic -Wall -O2 -funroll-loops -mpopcnt -DFLINT_CPIMPORT=\"/userdata/raoulb/lib/share/flint/CPimport.txt\"


INCS=-I$(CURDIR) -I$(GMP_INCLUDE_DIR) -I$(MPFR_INCLUDE_DIR) -I$(FLINT_INCLUDE_DIR)

LIBS=-L$(CURDIR) -L$(FLINT_LIB_DIR) -L$(GMP_LIB_DIR) -L$(MPFR_LIB_DIR) -larb -lflint -lmpfr -lgmp -lpthread -lm

bkes:
	$(CC) $(CFLAGS) $(INCS) libkes2.h kes.c $(LIBS) -o kes

ekes:
	$(CC) $(CFLAGS) $(INCS) libkes2.h kes_enumerate.c $(LIBS) -o kes_enumerate

rekes:
	$(CC) $(CFLAGS) $(INCS) libkes2.h kes_rec_enumerate.c $(LIBS) -o kes_rec_enumerate
