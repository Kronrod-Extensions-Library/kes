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


void evaluate_weights_formula_legendre(acb_ptr, const acb_ptr, const int, long);
void evaluate_weights_formula_laguerre(acb_ptr, const acb_ptr, const int, long);
void evaluate_weights_formula_hermite_pro(acb_ptr, const acb_ptr, const int, long);


void evaluate_weights_formula_legendre(acb_ptr weights,
				       const acb_ptr nodes,
				       const int n,
				       long prec) {
    /* Compute the Gauss-Legendre quadrature weights by the analytic formula.
     *
     */
    int k;
    arb_t pf;

    fmpq_poly_t P;
    acb_ptr t;

    // The prefactor
    // 2 / n^2
    arb_init(pf);

    arb_set_si(pf, n);
    arb_pow_ui(pf, pf, 2, prec);
    arb_ui_div(pf, 2, pf, prec);

    // The other part
    // (gamma^2 - 1) / (gamma * P_n(gamma) - P_{n-1}(gamma))^2
    fmpq_poly_init(P);
    t = _acb_vec_init(n);

    legendre_polynomial(P, n);
    evaluate_polynomial_vector(t, P, nodes, n, prec);
    legendre_polynomial(P, n-1);
    evaluate_polynomial_vector(weights, P, nodes, n, prec);

    for(k = 0; k < n; k++) {
        acb_mul((t+k), (t+k), (nodes+k), prec);
	acb_sub((weights+k), (t+k), (weights+k), prec);
	acb_pow_ui((weights+k), (weights+k), 2, prec);
	acb_one((t+k));
	acb_submul((t+k), (nodes+k), (nodes+k), prec);
        acb_div((weights+k), (t+k), (weights+k), prec);
	acb_mul_arb((weights+k), (weights+k), pf, prec);
    }

    arb_clear(pf);
    _acb_vec_clear(t, n);
    fmpq_poly_clear(P);
}


void evaluate_weights_formula_laguerre(acb_ptr weights,
				       const acb_ptr nodes,
				       const int n,
				       long prec) {
    /* Compute the Gauss-Laguerre quadrature weights by the analytic formula.
     *
     */
    int k;
    arb_t pf;

    fmpq_poly_t L;

    // The prefactor
    // 1 / (n+1)^2
    arb_init(pf);

    arb_set_si(pf, n);
    arb_add_si(pf, pf, 1, prec);
    arb_pow_ui(pf, pf, 2, prec);
    arb_inv(pf, pf, prec);

    // The other part
    // gamma / L^2_{n+1}(gamma)
    fmpq_poly_init(L);

    laguerre_polynomial(L, n+1);
    evaluate_polynomial_vector(weights, L, nodes, n, prec);

    for(k = 0; k < n; k++) {
        acb_pow_ui((weights+k), (weights+k), 2, prec);
        acb_div((weights+k), (nodes+k), (weights+k), prec);
        acb_mul_arb((weights+k), (weights+k), pf, prec);
    }

    arb_clear(pf);
    fmpq_poly_clear(L);
}


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
