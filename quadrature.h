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

#include "numerics.h"


void evaluate_weights_formula_hermite(acb_ptr, const acb_ptr, const int, long);


void evaluate_weights_formula_hermite_pro(acb_ptr weights,
                                          const acb_ptr nodes,
                                          const int n,
                                          long prec) {
    /* Compute the Gauss-Hermite quadrature weights by the analytic formula.
     *
     */
    int k;
    arb_t t;
    arb_t pf;

    fmpq_poly_t H;

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
    fmpq_poly_init(H);

    hermite_polynomial_pro(H, n-1);
    evaluate_polynomial_vector(weights, H, nodes, n, prec);

    for(k = 0; k < n; k++) {
        acb_mul((weights+k), (weights+k), (weights+k), prec);
        acb_inv((weights+k), (weights+k), prec);
        acb_mul_arb((weights+k), (weights+k), pf, prec);
    }

    arb_clear(t);
    arb_clear(pf);
    fmpq_poly_clear(H);
}


#endif
