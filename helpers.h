/*  Author: R. Bourquin
 *  Copyright: (C) 2014 R. Bourquin
 *  License: GNU GPL v2 or above
 *
 *  A library of helper functions to search for
 *  Kronrod extensions of Gauss quadrature rules.
 */

#ifndef __HH__helpers
#define __HH__helpers

void ps(const int, const int, const int);
void logit(const int, const int, const char *, ...);


void ps(const int loglevel, const int verbosity, const int n) {
    int i;
    if(verbosity >= loglevel) {
        for(i = 0; i < n; i++) {
            printf("  ");
        }
    }
}


void logit(const int loglevel, const int verbosity, const char *format, ...) {
    #if PRINTLOG
    va_list argp;
    if(verbosity >= loglevel) {
        va_start(argp, format);
        vprintf(format, argp);
        va_end(argp);
    }
    #endif
}

#endif
