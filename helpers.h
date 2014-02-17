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
int logit(const int, const int, const char *, ...);


void ps(const int loglevel, const int verbosity, const int n) {
    int i;
    if(verbosity >= loglevel) {
        for(i = 0; i < n; i++) {
            printf("  ");
        }
    }
}


int logit(const int loglevel, const int verbosity, const char *format, ...) {
    va_list argp;
    int ret;
    ret = 0;
    if(verbosity >= loglevel) {
        va_start(argp, format);
        ret = vprintf(format, argp);
        va_end(argp);
    }
    return ret;
}

#endif
