/*  Author: R. Bourquin
 *  Copyright: (C) 2014 R. Bourquin
 *  License: GNU GPL v2 or above
 *
 *  A library of helper functions to search for
 *  Kronrod extensions of Gauss quadrature rules.
 */

#ifndef __HH__quadrature
#define __HH__quadrature

#include <stdarg.h>

#include "arb.h"
#include "acb.h"


#define NCHECKDIGITS 53


void evaluate_hermite(acb_t, const acb_t, const int, const long);
void evaluate_hermite_vector(acb_ptr, const acb_ptr, const int, const int, long);
void evaluate_weights_formula_hermite(acb_ptr, const acb_ptr, const int, long);


void evaluate_hermite(acb_t value,
                      const acb_t node,
                      const int n,
                      const long prec) {
    /*
     *
     */
    int i;

    acb_t H0;
    acb_t H1;
    acb_t x;
    acb_t T0;
    acb_t T1;

    acb_init(x);
    acb_init(H0);
    acb_init(H1);

    /* This is x */
    acb_mul_si(x, node, 1, prec);

    acb_one(H0);
    acb_set(H1, x);

    if(n == 0) {
        acb_set(value, H0);
    } else if(n == 1) {
        acb_set(value, H1);
    } else {
        acb_init(T0);
        acb_init(T1);
        for(i = 1; i < n; i++) {
            acb_mul(T1, x, H1, prec);
            acb_mul_si(T0, H0, i, prec);
            acb_sub(value, T1, T0, prec);
            acb_set(H0, H1);
            acb_set(H1, value);
        }
        acb_clear(T0);
        acb_clear(T1);
    }

    acb_clear(H0);
    acb_clear(H1);
    acb_clear(x);
}


void evaluate_hermite_vector(acb_ptr values,
                             const acb_ptr nodes,
                             const int N,
                             const int K,
                             long prec) {
    /*
     *
     */
    int k;
    for(k = 0; k < K; k++) {
        evaluate_hermite((values+k), nodes+k, N, prec);
    }
}

void evaluate_weights_formula_hermite(acb_ptr weights,
                                      const acb_ptr nodes,
                                      const int n,
                                      long prec) {
    /*
     *
     */
    int k;
    arb_t t;
    arb_t pf;

    // The prefactor
    // 2^(n-1) Gamma(n+1) sqrt(pi) / n^2  prob version
    // 2^(-1) Gamma(n+1) / n^2            phys version
    arb_init(t);
    arb_init(pf);

    arb_ui_pow_ui(pf, 2, 0, prec);
    arb_fac_ui(t, n, prec);
    arb_mul(pf, pf, t, prec);
    arb_div_ui(pf, pf, n*n, prec);

    // The other part
    // 1 / H^2_{n-1}(gamma)
    evaluate_hermite_vector(weights, nodes, n-1, n, prec);

    for(k = 0; k < n; k++) {
        acb_mul((weights+k), (weights+k), (weights+k), prec);
        acb_inv((weights+k), (weights+k), prec);
        acb_mul_arb((weights+k), (weights+k), pf, prec);
    }

    arb_clear(t);
    arb_clear(pf);
}


#endif
