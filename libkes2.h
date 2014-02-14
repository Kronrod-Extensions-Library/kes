/*  Author: R. Bourquin
 *  Copyright: (C) 2014 R. Bourquin
 *  License: GNU GPL v2 or above
 *
 *  A library of helper functions to search for
 *  Kronrod extensions of Gauss quadrature rules.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpq_poly.h"
#include "flint/fmpz_mat.h"
#include "flint/fmpq_mat.h"

#include "fmprb.h"
#include "fmpcb_poly.h"
#include "fmpcb_mat.h"
/*#include "arith.h"*/


#define NCHECKDIGITS 20


int find_extension(fmpq_poly_t En, const fmpq_poly_t Pn, const int nu, const int loglevel);
int find_multi_extension(fmpq_poly_t, const fmpq_poly_t, const int, const int[], const int, const int);

void recursive_enumerate(const fmpq_poly_t, const int, const int, const int, const int, fmpz_mat_t, const int);

inline void compute_nodes(fmpcb_ptr, const fmpq_poly_t, const long, const int);
void compute_nodes_and_weights(fmpcb_ptr, fmpcb_ptr, const fmpq_poly_t, const long, const int);

int validate_rule(long*, long*, const fmpq_poly_t, const long, const int);
int validate_extension_by_poly(long*, const fmpq_poly_t, const long, const int);
int validate_extension_by_roots(const fmpcb_ptr, const long, const long, const int);
int validate_extension_by_weights(const fmpcb_ptr, const long, const long, const int);
long validate_real_roots(const fmpcb_ptr, const long, const long, const int);
long validate_positive_weights(const fmpcb_ptr, const long, const long, const int);

void integrate(fmpq_t M, const int n);
void integrate_hermite_pro(fmpq_t I, const int n);
void integrate_hermite_phy(fmpq_t I, const int n);
void hermite_polynomial_pro(fmpq_poly_t Hn, const int n);
void hermite_polynomial_phy(fmpq_poly_t Hn, const int n);
/*
void laguerre_polynomial(fmpq_poly_t Ln, const int n);
void legendre_polynomial(fmpq_poly_t Ln, const int n);
*/

void poly_roots(fmpcb_ptr, const fmpq_poly_t, const long, const long, const int);
int check_accuracy(const fmpcb_ptr, const long, const long);

void ps(const int, const int, const int);
int logit(const int, const int, const char *, ...);


/**********************************************************************/


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


int find_extension(fmpq_poly_t En,
                   const fmpq_poly_t Pn,
                   const int p,
                   const int loglevel) {
    /* Extend the degree n polynomial Pk by one Kronrod extension Ek of degree p.
     *
     * We search for a monic polynomial G(x) such that
     * \int_\Omega F(t) G(t) t^i \rho(t) dt = 0   for all   i = 0, ..., p-1.
     * To obtain the coefficients  a_0, ..., a_{p-1} of G(t)  we try to solve
     * a  p x p  linear system. If successful the extension E_k exists.
     *
     * Note that the extension E_n is unly valid if E_k has p real roots
     * inside the domain \Omega and if further all weights are positive.
     *
     * En: The polynomial defining the extension
     * Pn: The polynomial defining the basis
     * p: The degree of the extension
     * loglevel: The log verbosity
     */
    slong deg_P;
    slong rows;
    slong cols;

    fmpq_mat_t M, rhs, X;
    int i, k, j;

    slong deg;
    fmpq_t coeff;
    fmpq_t integral;
    fmpq_t element;

    int solvable;

    deg_P = fmpq_poly_degree(Pn);

    fmpq_init(coeff);
    fmpq_init(integral);
    fmpq_init(element);

    rows = p;
    cols = p;
    fmpq_mat_init(M, rows, cols);
    fmpq_mat_zero(M);
    fmpq_mat_init(rhs, rows, 1);
    fmpq_mat_zero(rhs);

    /* Build the system matrix  */
    for(i = 0; i < rows; i++) {
        for(k = 0; k < cols; k++) {
            fmpq_zero(element);
            for(j = 0; j <= deg_P; j++) {
                fmpq_poly_get_coeff_fmpq(coeff, Pn, j);
                deg = i + k + j;
                integrate(integral, deg);
                fmpq_mul(integral, coeff, integral);
                fmpq_add(element, element, integral);
            }
            fmpq_set(fmpq_mat_entry(M, i, k), element);
        }
    }

    /* Build the right hand side */
    for(i = 0; i < rows; i++) {
        fmpq_zero(element);
        for(j = 0; j <= deg_P; j++) {
            fmpq_poly_get_coeff_fmpq(coeff, Pn, j);
            deg = i + p + j;
            integrate(integral, deg);
            fmpq_mul(integral, coeff, integral);
            fmpq_add(element, element, integral);
        }
        fmpq_set(fmpq_mat_entry(rhs, i, 0), element);
    }
    fmpq_mat_neg(rhs, rhs);

    /*fmpq_mat_print(M);
      fmpq_mat_print(rhs);*/

    /* Try to solve the system */
    fmpq_mat_init(X, rows, 1);
    fmpq_mat_zero(X);
    solvable = fmpq_mat_solve_fraction_free(X, M, rhs);

    /*fmpq_mat_print(X);*/
    logit(1, loglevel, "Solvable: %i\n", solvable);

    /* Assemble the polynomial */
    fmpq_poly_zero(En);
    if(solvable) {
        for(k = 0; k < rows; k++) {
            fmpq_poly_set_coeff_fmpq(En, k, fmpq_mat_entry(X, k, 0));
        }
        fmpq_poly_set_coeff_si(En, p, 1);
    }

    fmpq_poly_canonicalise(En);

    /* Clean up */
    fmpq_clear(coeff);
    fmpq_clear(integral);
    fmpq_clear(element);
    fmpq_mat_clear(M);
    fmpq_mat_clear(rhs);
    fmpq_mat_clear(X);
    return solvable;
}


int find_multi_extension(fmpq_poly_t En,
                         const fmpq_poly_t P,
                         const int k,
                         const int levels[],
                         const int validate_extension,
                         const int loglevel) {
    /* Extends the polynomial P repeatedly by nested Kronrod extensions.
     * En: The polynomial defining the extension
     * P: The polynomial defining the basic rule
     * k: The number of nested extensions
     * levels: An array with the extension levels p_1, ..., p_k
     * validate_extension: Check if the roots are all real
     * loglevel: The log verbosity
     */
    int i;
    int solvable, valid;
    fmpq_poly_t Pn;
    long nrroots;
    char *strf;
    int success;

    fmpq_poly_init(Pn);
    fmpq_poly_set(Pn, P);
    strf = fmpq_poly_get_str_pretty(P, "t");

    success = 1;

    for(i = 1; i < k; i++) {
        logit(1, loglevel, "-------------------------------------------------\n");
        logit(1, loglevel, "Trying to find an order %i Kronrod extension for:\n", levels[i]);
        if(loglevel >= 2) {
            strf = fmpq_poly_get_str_pretty(Pn, "t");
            flint_printf("P%i : %s\n", i, strf);
        }

        solvable = find_extension(En, Pn, levels[i], loglevel);

        if(!solvable) {
            success = 0;
            printf("******************************\n");
            printf("*** EXTENSION NOT SOVLABLE ***\n");
            printf("******************************\n");
            break;
        }

        if(validate_extension) {
            valid = validate_extension_by_poly(&nrroots, En, NCHECKDIGITS, loglevel);

            if(!valid) {
                success = 0;
                printf("************************************\n");
                printf("*** EXTENSION WITH INVALID NODES ***\n");
                printf("************************************\n");
                break;
            }
        }

        /* Iterate */
        fmpq_poly_mul(Pn, Pn, En);
        fmpq_poly_canonicalise(Pn);
    }

    flint_free(strf);
    fmpq_poly_clear(Pn);
    return success;
}


void recursive_enumerate(const fmpq_poly_t Pn,
                         const int n,
                         const int maxp,
                         const int rec,
                         const int maxrec,
                         fmpz_mat_t table,
                         const int loglevel) {
    /* Recursively enumerate quadrature rules
     */
    int p;
    int solvable, valid;
    long nrroots;
    fmpq_poly_t Pnp1, Enp1, En;

    ps(1, loglevel, rec);
    logit(1, loglevel, "Trying to find extension of (on layer %i):\n", rec);

    fmpq_poly_init(Pnp1);
    fmpq_poly_init(En);
    fmpq_poly_init(Enp1);

    /* Loop over possible (non-recursive) extensions  */
    for(p = 1/*(n-1)/2+1*/; p <= maxp; p++) {

        solvable = find_extension(En, Pn, p, loglevel);
        valid = validate_extension_by_poly(&nrroots, En, NCHECKDIGITS, loglevel);

        if(solvable && valid) {
            ps(1, loglevel, rec);
            logit(1, loglevel, "Found valid extension for n: %i and p: %i (on layer %i)\n", n, p, rec);

            printf("RULE: %i  ", rec+2);
            fmpz_set_ui(fmpz_mat_entry(table, rec+1, 0), p);
            fmpz_mat_print(table);
            printf("\n");
            /*stre = fmpq_poly_get_str_pretty(En, "t");
            ps(rec);
            flint_printf("E%i : %s\n", n, stre);*/

            /* Follow the recursion down */
            if(rec <= maxrec) {
                ps(1, loglevel, rec);
                logit(1, loglevel, "==> Going down, new layer: %i\n", rec+1);
                fmpq_poly_mul(Pnp1, Pn, En);
                recursive_enumerate(Pnp1, n+p, maxp, rec+1, maxrec, table, loglevel);
            } else {
                ps(1, loglevel, rec);
                logit(1, loglevel, "##> Maximum recursion depth reached, not descending\n");
            }

        } else {
            ps(1, loglevel, rec);
            logit(1, loglevel, "No valid extension for n: %i and p: %i found (on layer %i)\n", n, p, rec);
        }
    }
    ps(1, loglevel, rec);
    logit(1, loglevel, "Maximal extension order p: %i reached\n", maxp);

    ps(1, loglevel, rec);
    logit(1, loglevel, "==> Going up, leaving layer: %i\n", rec);

    fmpz_set_ui(fmpz_mat_entry(table, rec+1, 0), 0);

    fmpq_poly_clear(Pnp1);
    fmpq_poly_clear(En);
    fmpq_poly_clear(Enp1);
    return;
}


inline void compute_nodes(fmpcb_ptr nodes,
                   const fmpq_poly_t poly,
                   const long prec,
                   const int loglevel) {
    /*
     * nodes: An array containing the nodes
     * poly: The polynomial whose roots define the nodes
     * prec: Number of bits in target precision
     * loglevel: The log verbosity
     */
    logit(1, loglevel, "-------------------------------------------------\n");
    logit(1, loglevel, "Computing nodes\n");
    poly_roots(nodes, poly, 53, prec, loglevel);
}


void compute_nodes_and_weights(fmpcb_ptr nodes,
                               fmpcb_ptr weights,
                               const fmpq_poly_t poly,
                               const long target_prec,
                               const int loglevel) {
    /*
     * nodes: An array containing the nodes
     * weights: An array containing the weights
     * poly: The polynomial whose roots define the nodes
     * target_prec: Number of bits in target precision
     * loglevel: The log verbosity
     */
    int k, j;
    slong K;
    fmpq_t integral;
    fmpcb_t element;
    fmpcb_mat_t A;
    fmpcb_mat_t B;
    fmpcb_mat_t X;
    int solvable;
    long initial_prec, prec;

    K = fmpq_poly_degree(poly);
    fmpcb_mat_init(A, K, K);
    fmpcb_mat_init(B, K, 1);
    fmpcb_mat_init(X, K, 1);
    fmpcb_mat_zero(A);
    fmpcb_mat_zero(B);
    fmpcb_mat_zero(X);

    logit(1, loglevel, "-------------------------------------------------\n");
    logit(1, loglevel, "Computing nodes and weights\n");

    /* Precision in number of bits */
    initial_prec = 53;

    for(prec = initial_prec; ; prec *= 2) {
        /* Find the roots up to prec bits */
        poly_roots(nodes, poly, prec, prec, loglevel);

        /* Build the system matrix */
        for(k = 0; k < K; k++) {
            for(j = 0; j < K; j++) {
                fmpcb_pow_ui(element, nodes+j, k, prec);
                fmpcb_set(fmpcb_mat_entry(A, k, j), element);
            }
        }

        /* Build the right hand side */
        for(k = 0; k < K; k++) {
            integrate(integral, k);
            fmpcb_set_fmpq(element, integral, prec);
            fmpcb_set(fmpcb_mat_entry(B, k, 0), element);
        }

        /* Solve the system and obtain weights */
        logit(4, loglevel, " current precision for weights: %ld\n", prec);

        solvable = fmpcb_mat_solve(X, A, B, prec);

        for(k = 0; k < K; k++) {
            *(k + weights) = *fmpcb_mat_entry(X, k, 0);
        }

        logit(4, loglevel, "Linear system for weights solvable: %i\n", solvable);

        /* Check accuracy of weights here */
        if(solvable && check_accuracy(weights, K, target_prec)) {
            logit(4, loglevel, "Sufficient bits for target precision reached\n");
            break;
        }
    }
}


int validate_rule(long* nrroots,
		  long* nnnweights,
		  const fmpq_poly_t En,
		  const long prec,
		  const int loglevel) {
    /*
     * nrroots: Number of real roots found
     * nnnweights: Number of non-negative weights found
     * En: Polynomial defining the extension
     * prec: Number of bits in target precision
     * loglevel: The log verbosity
     */
    slong deg;
    fmpcb_ptr roots, weights;
    long rroots, nnweights;

    /* This extension is invalid */
    deg = fmpq_poly_degree(En);
    if(deg < 0) {
        return 0;
    }

    /* Compute the roots */
    roots = _fmpcb_vec_init(deg);
    weights = _fmpcb_vec_init(deg);

    compute_nodes_and_weights(roots, weights, En, prec, loglevel);

    /* Validate roots */
    rroots = validate_real_roots(roots, deg, prec, loglevel);
    (*nrroots) = rroots;

    /* Validate weights */
    nnweights = validate_positive_weights(weights, deg, prec, loglevel);
    (*nnnweights) = nnweights;

    _fmpcb_vec_clear(roots, deg);
    _fmpcb_vec_clear(weights, deg);

    return (rroots == deg && nnweights == deg) ? 1 : 0;
}


int validate_extension_by_poly(long* nrroots,
                               const fmpq_poly_t En,
                               const long prec,
                               const int loglevel) {
    /*
     * nrroots: Number of real roots found
     * En: Polynomial defining the extension
     * prec: Number of bits in target precision
     * loglevel: The log verbosity
     */
    slong deg;
    fmpcb_ptr roots;
    long valid_roots;

    /* This extension is invalid */
    deg = fmpq_poly_degree(En);
    if(deg <= 0) {
        return 0;
    }

    /* Compute the roots */
    roots = _fmpcb_vec_init(deg);
    compute_nodes(roots, En, prec, loglevel);

    /* Validate roots */
    valid_roots = validate_real_roots(roots, deg, prec, loglevel);
    (*nrroots) = valid_roots;

    logit(1, loglevel, "Extension rule has valid nodes: %i\n", valid_roots == deg ? 1 : 0);

    _fmpcb_vec_clear(roots, deg);

    return valid_roots == deg ? 1 : 0;
}


int validate_extension_by_roots(const fmpcb_ptr roots,
                                const long deg,
                                const long prec,
                                const int loglevel) {
    /*
     * roots: Array of all roots of the extension
     * deg: Number of roots in the input array
     * prec: Number of bits in target precision
     * loglevel: The log verbosity
     */
    long valid_roots;

    /* This extension is invalid */
    if(deg <= 0) {
        return 0;
    }

    /* Validate roots */
    valid_roots = validate_real_roots(roots, deg, prec, loglevel);

    logit(1, loglevel, "Extension rule has valid nodes: %i\n", valid_roots == deg ? 1 : 0);

    return valid_roots == deg ? 1 : 0;
}


int validate_extension_by_weights(const fmpcb_ptr weights,
				  const long deg,
				  const long prec,
				  const int loglevel) {
    /*
     * weights: Array of all weights of the extension
     * deg: Number of weights in the input array
     * prec: Number of bits in target precision
     * loglevel: The log verbosity
     */
    long valid_weights;

    /* This extension is invalid */
    if(deg <= 0) {
        return 0;
    }

    /* Validate weights */
    valid_weights = validate_positive_weights(weights, deg, prec, loglevel);

    logit(1, loglevel, "Extension rule has valid weights: %i\n", valid_weights == deg ? 1 : 0);

    return valid_weights == deg ? 1 : 0;
}


long validate_real_roots(const fmpcb_ptr roots,
                         const long n,
                         const long prec,
                         const int loglevel) {
    /* Validate the roots found
     * This part is heuristic: we check if the imaginary part is small enough,
     *                         but this does not guarantee the root is real.
     *
     * roots: Array containing the roots to validate
     * n: Number of roots in the input array
     * prec: The number of bits used for validation.
     * loglevel: The log verbosity
     */
    long valid_roots;
    long i;

    valid_roots = 0;
    for(i = 0; i < n; i++) {
        if(fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_imagref(roots + i)), -prec) <= 0) {
            valid_roots++;
        }
    }

    logit(1, loglevel, "Valid roots found: %lu out of %lu\n", valid_roots, n);

    return valid_roots;
}


long validate_positive_weights(const fmpcb_ptr weights,
                               const long n,
                               const long prec,
                               const int loglevel) {
    /* Validate the weights found
     * This part is heuristic: we check if the imaginary part is small enough,
     *                         but this does not guarantee the weight is real.
     *                         Additionally we check that the ball is not negative,
     *                         but this does not guarantee the weight is positive.
     *
     * weights: Array containing the weights to validate
     * n: Number of weights in the input array
     * prec: The number of bits used for validation.
     * loglevel: The log verbosity
     */
    long positive_weights;
    long negative_weights;
    long indefinite_weights;
    long i;

    positive_weights = 0;
    negative_weights = 0;
    indefinite_weights = 0;
    for(i = 0; i < n; i++) {
        if(fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_imagref(weights + i)), -prec) <= 0) {
            if(fmprb_is_positive(fmpcb_realref(weights + i))) {
                positive_weights++;
            } else if(fmprb_is_negative(fmpcb_realref(weights + i))) {
                negative_weights++;
            } else {
                indefinite_weights++;
            }
        }
    }

    logit(1, loglevel, "Positive weights:   %ld out of %ld\n", positive_weights, n);
    logit(1, loglevel, "Indefinite weights: %ld out of %ld\n", indefinite_weights, n);
    logit(1, loglevel, "Negative weights:   %ld out of %ld\n", negative_weights, n);

    return positive_weights + indefinite_weights;
}








void integrate(fmpq_t M, const int n) {
    integrate_hermite_pro(M, n);
}


void integrate_hermite_pro(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{-\infty}^\infty \exp(-x^2/2) x^n dx
     *
     * n even:  I = \Gamma{\frac{n+1}{2}}
     * n odd:   I = 0
     *
     * 1  0  1  0  3  0  15  0  105  0  945  0  10395  0  135135
     * We omit a factor of \sqrt{2\pi}
     */
    int i;
    fmpz_t tmp;

    fmpz_init(tmp);

    if(n % 2 == 1) {
        fmpq_zero(I);
    } else {
        fmpq_one(I);
        for(i = 1; i <= n/2; i++) {
            fmpz_set_ui(tmp, 2*i - 1);
            fmpq_mul_fmpz(I, I, tmp);
        }
    }

    fmpz_clear(tmp);
}


void integrate_hermite_phy(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{-\infty}^\infty \exp(-x^2/2) x^n dx
     *
     * n even:  I = 2^{\frac{n+1}{2}} \Gamma{\frac{n+1}{2}}
     * n odd:   I = 0
     *
     * 1  0  1/2  0  3/4  0  15/8  0  105/16  0  945/32  0  10395/64  0  135135/128
     * We omit a factor of \sqrt{\pi}
     */
    int i;
    fmpz_t tmp;

    fmpz_init(tmp);

    if(n % 2 == 1) {
        fmpq_zero(I);
    } else {
        fmpq_one(I);
        for(i = 1; i <= n/2; i++) {
            fmpz_set_ui(tmp, 2*i - 1);
            fmpq_mul_fmpz(I, I, tmp);
        }
        fmpq_div_2exp(I, I, n/2);
    }

    fmpz_clear(tmp);
}


void hermite_polynomial_pro(fmpq_poly_t Hn, const int n) {
    /* Compute the n-th Hermite polynomial by a
     * three term recursion:
     *
     * H_0(x) = 1
     * H_1(x) = x
     *
     * H_{n+1}(x) = x H_n(x) - n H_{n-1}(x)
     *
     * This sequence yields the probabilists' Hermite polynomials.
     */
    fmpq_poly_t H0, H1;
    fmpq_poly_t x;
    fmpq_poly_t T0, T1;
    int i;

    fmpq_poly_init(H0);
    fmpq_poly_init(H1);
    fmpq_poly_init(x);

    fmpq_poly_set_str(x, "2  0 1");
    fmpq_poly_canonicalise(x);

    fmpq_poly_set_str(H0, "1  1");
    fmpq_poly_canonicalise(H0);

    fmpq_poly_set(H1, x);
    fmpq_poly_canonicalise(H1);

    if(n == 0) {
        fmpq_poly_set(Hn, H0);
    } else if(n == 1) {
        fmpq_poly_set(Hn, H1);
    } else {
        fmpq_poly_init(T0);
        fmpq_poly_init(T1);
        for(i = 1; i < n; i++) {
            fmpq_poly_mul(T1, x, H1);
            fmpq_poly_scalar_mul_si(T0, H0, i);
            fmpq_poly_sub(Hn, T1, T0);
            fmpq_poly_set(H0, H1);
            fmpq_poly_set(H1, Hn);
        }
        fmpq_poly_clear(T0);
        fmpq_poly_clear(T1);
    }

    fmpq_poly_clear(H0);
    fmpq_poly_clear(H1);

    fmpq_poly_canonicalise(Hn);
    return;
}


void hermite_polynomial_phy(fmpq_poly_t Hn, const int n) {
    /* Compute the n-th Hermite polynomial by a
     * three term recursion:
     *
     * H_0(x) = 1
     * H_1(x) = x
     *
     * H_{n+1}(x) = 2 x H_n(x) - 2 n H_{n-1}(x)
     *
     * This sequence yields the physicists' Hermite polynomials.
     */
    fmpq_poly_t H0, H1;
    fmpq_poly_t x;
    fmpq_poly_t T0, T1;
    int i;

    fmpq_poly_init(H0);
    fmpq_poly_init(H1);
    fmpq_poly_init(x);

    /* This is 2x */
    fmpq_poly_set_str(x, "2  0 2");
    fmpq_poly_canonicalise(x);

    fmpq_poly_set_str(H0, "1  1");
    fmpq_poly_canonicalise(H0);

    fmpq_poly_set(H1, x);

    if(n == 0) {
        fmpq_poly_set(Hn, H0);
    } else if(n == 1) {
        fmpq_poly_set(Hn, H1);
    } else {
        fmpq_poly_init(T0);
        fmpq_poly_init(T1);
        for(i = 1; i < n; i++) {
            fmpq_poly_mul(T1, x, H1);
            fmpq_poly_scalar_mul_si(T0, H0, 2*i);
            fmpq_poly_sub(Hn, T1, T0);
            fmpq_poly_set(H0, H1);
            fmpq_poly_set(H1, Hn);
        }
        fmpq_poly_clear(T0);
        fmpq_poly_clear(T1);
    }

    fmpq_poly_clear(H0);
    fmpq_poly_clear(H1);

    return;
}


/* The following function are borrowed from the "poly_roots"
   example in the arb library documentation. The header of
   the "poly_roots.c" files states:

       This file is public domain. Author: Fredrik Johansson.
*/


void poly_roots(fmpcb_ptr roots,
                const fmpq_poly_t poly,
                const long initial_prec,
                const long target_prec,
                const int loglevel) {
    /*
     * roots: An array containing the roots
     * poly: The polynomial whose roots to compute
     * initial_prec: Number of bits in initial precision
     * target_prec: Number of bits in target precision
     * loglevel: The log verbosity
     */
    long prec, deg, isolated, maxiter;
    fmpcb_poly_t cpoly;

    deg = fmpq_poly_degree(poly);
    fmpcb_poly_init(cpoly);

    for(prec = initial_prec; ; prec *= 2) {
        fmpcb_poly_set_fmpq_poly(cpoly, poly, prec);
        maxiter = FLINT_MIN(deg, prec);

        logit(4, loglevel, "  current precision for roots: %ld\n", prec);
        isolated = fmpcb_poly_find_roots(roots, cpoly, prec == initial_prec ? NULL : roots, maxiter, prec);

        if(isolated == deg && check_accuracy(roots, deg, target_prec)) {
            break;
        }
    }
    fmpcb_poly_clear(cpoly);
}


int check_accuracy(const fmpcb_ptr vec, const long len, const long prec) {
    /* Check if all balls in a vector have a radius small enough
     * to fit the target precision.
     *
     * vec: Vector of balls to test
     * len: Number of balls to test
     * prec: Target precision in number bits
     */
    long i;

    for(i = 0; i < len; i++) {
        if(   fmpr_cmp_2exp_si(fmprb_radref(fmpcb_realref(vec + i)), -prec) >= 0
           || fmpr_cmp_2exp_si(fmprb_radref(fmpcb_imagref(vec + i)), -prec) >= 0)
            return 0;
    }
    return 1;
}
