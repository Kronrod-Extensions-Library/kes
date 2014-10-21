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


void sort_nodes(acb_ptr, const int);
void evaluate_weights_formula_legendre(acb_ptr, const acb_ptr, const int, const long);
void evaluate_weights_formula_laguerre(acb_ptr, const acb_ptr, const int, const long);
void evaluate_weights_formula_hermite_pro(acb_ptr, const acb_ptr, const int, const long);
void evaluate_weights_formula_chebyshevt(acb_ptr, const acb_ptr, const int, const long);
void evaluate_weights_formula_chebyshevu(acb_ptr, const acb_ptr, const int, const long);


void sort_nodes(acb_ptr nodes, const int n) {
    /* Sort quadrature nodes in-place
     *
     * n: Number of nodes
     */
    int i, j;
    acb_t t;
    acb_init(t);
    // In-place insertion sort (for small n)
    for(i = 1; i < n; i++) {
        acb_set(t, (nodes+i));
        j = i;
        while(j > 0 &&
              arf_cmp(arb_midref(acb_realref( (nodes+j-1) )),
                      arb_midref(acb_realref( t ))) > 0) {
            acb_set((nodes+j), (nodes+j-1));
            j--;
        }
        acb_set((nodes+j), t);
    }
    acb_clear(t);
}


void evaluate_weights_formula_legendre(acb_ptr weights,
                                       const acb_ptr nodes,
                                       const int n,
                                       const long prec) {
    /* Compute the Gauss-Legendre quadrature weights by the analytic formula.
     *
     * w_k = \frac{2}{(1-x_k^2) P'_n(x_k)^2}
     * k = 0, ..., n-1
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
                                       const long prec) {
    /* Compute the Gauss-Laguerre quadrature weights by the analytic formula.
     *
     * w_k = \frac{x_k}{(n+1)^2 L_{n+1}(x_k)^2}
     * k = 0, ..., n-1
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
                                          const long prec) {
    /* Compute the Gauss-Hermite quadrature weights by the analytic formula.
     *
     * w_k = \frac{2^{n-1} n! \sqrt{\pi}}{n^2 H_{n-1}(x_k)^2}  phys version
     * w_k = \frac{n!}{n^2 H_{n-1}(x_k)^2}                     prob version
     * k = 0, ..., n-1
     */
    int k;
    arb_t pf;
    fmpq_poly_t H;

    // The prefactor
    // 2^(n-1) Gamma(n+1) sqrt(pi) / n^2  phys version
    // Gamma(n+1) / n^2                   prob version
    arb_init(pf);

    arb_fac_ui(pf, n, prec);
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

    arb_clear(pf);
    fmpq_poly_clear(H);
}


void evaluate_weights_formula_chebyshevt(acb_ptr weights,
                                         const acb_ptr nodes,
                                         const int n,
                                         const long prec) {
    /* Compute the Gauss-Chebyshev quadrature weights by the analytic formula.
     *
     * w_k = \frac{\pi}{n}
     * k = 0, ..., n-1
     */
    int k;

    for(k = 0; k < n; k++) {
        acb_const_pi((weights+k), prec);
        acb_div_ui((weights+k), (weights+k), n, prec);
    }
}


void evaluate_weights_formula_chebyshevu(acb_ptr weights,
                                         const acb_ptr nodes,
                                         const int n,
                                         const long prec) {
    /* Compute the Gauss-Chebyshev quadrature weights by the analytic formula.
     *
     * w_k = \frac{\pi}{n + 1} \sin( \frac{k+1}{n+1}\pi )^2
     * k = 0, ..., n-1
     */
    int k;
    arb_t t;
    acb_t pf;

    // The prefactor
    // pi / (n+1)
    acb_init(pf);

    acb_const_pi(pf, prec);
    acb_div_ui(pf, pf, n+1, prec);

    // The other part
    // sin( ((k+1)*pi) / (n+1) )^2
    arb_init(t);

    for(k = 0; k < n; k++) {
        arb_const_pi(t, prec);
        arb_mul_ui(t, t, k+1, prec);
        arb_div_ui(t, t, n+1, prec);
        arb_sin(t, t, prec);
        arb_pow_ui(t, t, 2, prec);
        acb_mul_arb((weights+k), pf, t, prec);
    }

    arb_clear(t);
    acb_clear(pf);
}


#endif
