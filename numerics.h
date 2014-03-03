/*  Author: R. Bourquin
 *  Copyright: (C) 2014 R. Bourquin
 *  License: GNU GPL v2 or above
 *
 *  A library of helper functions to search for
 *  Kronrod extensions of Gauss quadrature rules.
 */

#ifndef __HH__numerics
#define __HH__numerics

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

#include "helpers.h"


long validate_real_roots(const fmpcb_ptr, const long, const long, const int);
long validate_real_nonnegative_roots(const fmpcb_ptr, const long, const long, const int);
long validate_real_interval_roots(const fmpcb_ptr, const long, const long, const int);

long validate_positive_weights(const fmpcb_ptr, const long, const long, const int);

void poly_roots(fmpcb_ptr, const fmpq_poly_t, const long, const long, const int);
int check_accuracy(const fmpcb_ptr, const long, const long);


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
        if(! (fmprb_is_positive(fmpcb_imagref(roots+i)) || fmprb_is_negative(fmpcb_imagref(roots+i)))) {
            valid_roots++;
        }
    }

    logit(1, loglevel, "Valid roots found: %lu out of %lu\n", valid_roots, n);

    return valid_roots;
}


long validate_real_nonnegative_roots(const fmpcb_ptr roots,
                                     const long n,
                                     const long prec,
                                     const int loglevel) {
    /* Validate the roots found
     * This part is heuristic: we check if the imaginary part is small enough
     *                         and if the real part if non-negative.
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
        if(! (fmprb_is_positive(fmpcb_imagref(roots+i)) || fmprb_is_negative(fmpcb_imagref(roots+i)))
           && fmprb_is_nonnegative(fmpcb_realref(roots+i))) {
            valid_roots++;
        }
    }

    logit(1, loglevel, "Valid roots found: %lu out of %lu\n", valid_roots, n);

    return valid_roots;
}


long validate_real_interval_roots(const fmpcb_ptr roots,
                                  const long n,
                                  const long prec,
                                  const int loglevel) {
    /* Validate the roots found
     * This part is heuristic: we check if the imaginary part is small enough
     *                         and if the real part if non-negative.
     *                         but this does not guarantee the root is real.
     *
     * roots: Array containing the roots to validate
     * n: Number of roots in the input array
     * prec: The number of bits used for validation.
     * loglevel: The log verbosity
     */
    long valid_roots;
    long i;
    fmprb_t one;
    fmprb_t mone;
    fmprb_init(one);
    fmprb_init(mone);
    fmprb_one(one);
    fmprb_neg(mone, one);

    valid_roots = 0;
    for(i = 0; i < n; i++) {
        if(! (fmprb_is_positive(fmpcb_imagref(roots+i)) || fmprb_is_negative(fmpcb_imagref(roots+i)))) {
            if(    fmpr_cmpabs(fmprb_midref(fmpcb_realref(roots+i)), fmprb_midref(one)) <= 0
                || fmprb_contains(fmpcb_realref(roots+i), one)
                || fmprb_contains(fmpcb_realref(roots+i), mone)) {
                valid_roots++;
            }
        }
    }
    fmprb_clear(one);
    fmprb_clear(mone);

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
        if(fmpr_cmpabs_2exp_si(fmprb_midref(fmpcb_imagref(weights+i)), -prec) <= 0) {
            if(fmprb_is_positive(fmpcb_realref(weights+i))) {
                positive_weights++;
            } else if(fmprb_is_negative(fmpcb_realref(weights+i))) {
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
        if(   fmpr_cmp_2exp_si(fmprb_radref(fmpcb_realref(vec+i)), -prec) >= 0
           || fmpr_cmp_2exp_si(fmprb_radref(fmpcb_imagref(vec+i)), -prec) >= 0)
            return 0;
    }
    return 1;
}


#endif
