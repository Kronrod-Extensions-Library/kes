/*  Author: R. Bourquin
 *  Copyright: (C) 2014 R. Bourquin
 *  License: GNU GPL v2 or above
 *
 *  A library of helper functions to search for
 *  Kronrod extensions of Gauss quadrature rules.
 */

#ifndef __HH__polynomials
#define __HH__polynomials

#include <stdlib.h>

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpq_poly.h"
#include "flint/fmpz_mat.h"
#include "flint/fmpq_mat.h"


void hermite_polynomial_pro(fmpq_poly_t, const int);
void hermite_polynomial_phy(fmpq_poly_t, const int);
void integrate_hermite_pro(fmpq_t, const int);
void integrate_hermite_phy(fmpq_t, const int);
void moments_hermite_pro(fmpq_mat_t, const int);

void laguerre_polynomial(fmpq_poly_t, const int);
void integrate_laguerre(fmpq_t, const int);

void legendre_polynomial(fmpq_poly_t, const int);
void integrate_legendre(fmpq_t, const int);
void moments_legendre(fmpq_mat_t, const int);

void chebyshevt_polynomial(fmpq_poly_t Tn, const int n);
void integrate_chebyshevt(fmpq_t I, const int n);
void moments_chebyshevt(fmpq_mat_t, const int);

void chebyshevu_polynomial(fmpq_poly_t Un, const int n);
void integrate_chebyshevu(fmpq_t I, const int n);
void moments_chebyshevu(fmpq_mat_t, const int);


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

    fmpq_poly_canonicalise(Hn);
    return;
}


void integrate_hermite_pro(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{-\infty}^\infty \exp(-x^2/2) x^n dx
     *
     * n even:  I = \sqrt{2} 2^n \Gamma{\frac{n+1}{2}}
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
    return;
}


void integrate_hermite_phy(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{-\infty}^{\infty} \exp(-x^2) x^n dx
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
    return;
}


void moments_hermite_pro(fmpq_mat_t moments, const int n) {
    /*
     *
     */
    int i;
    fmpq_t tmp, entry;

    fmpq_mat_init(moments, 1, n);
    fmpq_init(tmp);
    fmpq_init(entry);
    fmpq_one(entry);

    for(i = 0; i < n; i++) {
        if(i % 2 == 0) {
            fmpq_set(fmpq_mat_entry(moments, 0, i), entry);
            fmpq_set_si(tmp, i + 1, 1);
            fmpq_mul(entry, entry, tmp);
        } else {
            fmpq_zero(fmpq_mat_entry(moments, 0, i));
        }
    }

    fmpq_clear(tmp);
    fmpq_clear(entry);
    return;
}


void laguerre_polynomial(fmpq_poly_t Ln, const int n) {
    /* Compute the n-th Laguerre polynomial by a
     * three term recursion:
     *
     * L_0(x) = 1
     * L_1(x) = 1 - x
     *
     * L_{n+1}(x) =   (2*n + 1 - x) / (n + 1) * L_n(x)
     *              - n / (n + 1) * L_{n-1}(x)
     *
     * This sequence yields Laguerre polynomials.
     */
    fmpq_poly_t L0, L1;
    fmpq_poly_t x;
    fmpq_poly_t T0, T1;
    int i;

    fmpq_poly_init(L0);
    fmpq_poly_init(L1);
    fmpq_poly_init(x);

    fmpq_poly_set_str(x, "2  1 -1");
    fmpq_poly_canonicalise(x);

    fmpq_poly_one(L0);
    fmpq_poly_canonicalise(L0);

    fmpq_poly_set(L1, x);
    fmpq_poly_canonicalise(L1);

    if(n == 0) {
        fmpq_poly_set(Ln, L0);
    } else if(n == 1) {
        fmpq_poly_set(Ln, L1);
    } else {
        fmpq_poly_init(T0);
        fmpq_poly_init(T1);
        for(i = 1; i < n; i++) {
            fmpq_poly_mul(T0, x, L1);
            fmpq_poly_scalar_mul_si(T1, L1, 2*i);
            fmpq_poly_add(T1, T1, T0);
            fmpq_poly_scalar_div_si(T1, T1, i+1);
            fmpq_poly_scalar_mul_si(T0, L0, i);
            fmpq_poly_scalar_div_si(T0, T0, i+1);
            fmpq_poly_sub(Ln, T1, T0);
            fmpq_poly_set(L0, L1);
            fmpq_poly_set(L1, Ln);
        }
        fmpq_poly_clear(T0);
        fmpq_poly_clear(T1);
    }

    fmpq_poly_clear(L0);
    fmpq_poly_clear(L1);

    fmpq_poly_canonicalise(Ln);
    return;
}


void integrate_laguerre(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{0}^\infty \exp(-x) x^n dx
     *
     * I = \Gamma{n+1}
     *
     * 1  1  2  6  24  120  720  5040  40320  362880  3628800
     */
    int i;
    fmpq_t t;

    fmpq_init(t);
    fmpq_one(I);
    for(i = 1; i <= n; i++) {
        fmpq_set_si(t, i, 1);
        fmpq_mul(I, I, t);
    }
    fmpq_clear(t);
    return;
}


void legendre_polynomial(fmpq_poly_t Pn, const int n) {
    /* Compute the n-th Legendre polynomial by a
     * three term recursion:
     *
     * P_0(x) = 1
     * P_1(x) = x
     *
     * P_{n+1}(x) = (2*n + 1) / (n + 1) * x * P_n(x)  -  n * P_{n-1}(x)
     *
     * This sequence yields Legendre polynomials.
     */
    fmpq_poly_t P0, P1;
    fmpq_poly_t x;
    fmpq_poly_t T0, T1;
    int i;

    fmpq_poly_init(P0);
    fmpq_poly_init(P1);
    fmpq_poly_init(x);

    fmpq_poly_set_str(x, "2  0 1");
    fmpq_poly_canonicalise(x);

    fmpq_poly_one(P0);
    fmpq_poly_canonicalise(P0);

    fmpq_poly_set(P1, x);
    fmpq_poly_canonicalise(P1);

    if(n == 0) {
        fmpq_poly_set(Pn, P0);
    } else if(n == 1) {
        fmpq_poly_set(Pn, P1);
    } else {
        fmpq_poly_init(T0);
        fmpq_poly_init(T1);
        for(i = 1; i < n; i++) {
            fmpq_poly_mul(T1, x, P1);
            fmpq_poly_scalar_mul_si(T1, T1, 2*i+1);
            fmpq_poly_scalar_div_si(T1, T1, i+1);
            fmpq_poly_scalar_mul_si(T0, P0, i);
            fmpq_poly_scalar_div_si(T0, T0, i+1);
            fmpq_poly_sub(Pn, T1, T0);
            fmpq_poly_set(P0, P1);
            fmpq_poly_set(P1, Pn);
        }
        fmpq_poly_clear(T0);
        fmpq_poly_clear(T1);
    }

    fmpq_poly_clear(P0);
    fmpq_poly_clear(P1);

    fmpq_poly_canonicalise(Pn);
    return;
}


void integrate_legendre(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{-1}^{1} x^n dx
     *
     * n even:  I = \frac{2}{n + 1}
     * n odd:   I = 0
     *
     * 2  0  2/3  0  2/5  0  2/7  0  2/9  0  2/11  0  2/13  0  2/15
     */
    fmpq_t t;

    fmpq_init(t);
    fmpq_one(I);

    if(n % 2 == 1) {
        fmpq_zero(I);
    } else {
        fmpq_set_si(t, 2, n+1);
        fmpq_mul(I, I, t);
    }
    fmpq_clear(t);
    return;
}


void moments_legendre(fmpq_mat_t moments, const int n) {
    /*
     *
     */
    int i;
    fmpq_t entry;

    fmpq_mat_init(moments, 1, n);
    fmpq_init(entry);
    fmpq_one(entry);

    for(i = 0; i < n; i++) {
        if(i % 2 == 0) {
            fmpq_set_si(entry, 2, i+1);
            fmpq_set(fmpq_mat_entry(moments, 0, i), entry);
        } else {
            fmpq_zero(fmpq_mat_entry(moments, 0, i));
        }
    }

    fmpq_clear(entry);
    return;
}


void chebyshevt_polynomial(fmpq_poly_t Tn, const int n) {
    /* Compute the n-th Chebyshev polynomial by a
     * three term recursion:
     *
     * T_0(x) = 1
     * T_1(x) = x
     *
     * T_{n+1}(x) = 2 * x * T_n(x)  -  T_{n-1}(x)
     *
     * This sequence yields Chebyshev polynomials.
     */
    fmpq_poly_t T0, T1;
    fmpq_poly_t x;
    int i;

    fmpq_poly_init(T0);
    fmpq_poly_init(T1);
    fmpq_poly_init(x);

    fmpq_poly_set_str(x, "2  0 1");
    fmpq_poly_canonicalise(x);

    fmpq_poly_one(T0);
    fmpq_poly_canonicalise(T0);

    fmpq_poly_set(T1, x);
    fmpq_poly_canonicalise(T1);

    if(n == 0) {
        fmpq_poly_set(Tn, T0);
    } else if(n == 1) {
        fmpq_poly_set(Tn, T1);
    } else {
        for(i = 1; i < n; i++) {
            fmpq_poly_mul(Tn, x, T1);
            fmpq_poly_scalar_mul_si(Tn, Tn, 2);
            fmpq_poly_sub(Tn, Tn, T0);
            fmpq_poly_set(T0, T1);
            fmpq_poly_set(T1, Tn);
        }
    }

    fmpq_poly_clear(T0);
    fmpq_poly_clear(T1);

    fmpq_poly_canonicalise(Tn);
    return;
}


void integrate_chebyshevt(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{-1}^{1} \frac{1}{\sqrt{1-x^2}} x^n dx
     *
     * n even:  I = \frac{2 \sqrt{\pi}}{n} \frac{\Gamma{\frac{n+1}{2}}}{\Gamma{\frac{n}{2}}}
     * n odd:   I = 0
     *
     * 1  0  1/2  0  3/8  0  5/16  0  35/128  0  63/256  0  231/1024  0  429/2048
     * We omit a factor of \pi
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
        for(i = 1; i <= n/2; i++) {
            fmpz_set_ui(tmp, i);
            fmpq_div_fmpz(I, I, tmp);
        }
        fmpq_div_2exp(I, I, n/2);
    }

    fmpz_clear(tmp);
    return;
}


void moments_chebyshevt(fmpq_mat_t moments, const int n) {
    /*
     *
     */
    int i, j;
    fmpz_t tmp;
    fmpq_t entry;

    fmpq_mat_init(moments, 1, n);
    fmpz_init(tmp);
    fmpq_init(entry);
    fmpq_one(entry);

    for(i = 0; i < n; i++) {
        if(i % 2 == 0) {
            fmpq_one(entry);
            for(j = 1; j <= i/2; j++) {
                fmpz_set_ui(tmp, 2*j - 1);
                fmpq_mul_fmpz(entry, entry, tmp);
            }
            for(j = 1; j <= i/2; j++) {
                fmpz_set_ui(tmp, j);
                fmpq_div_fmpz(entry, entry, tmp);
            }
            fmpq_div_2exp(fmpq_mat_entry(moments, 0, i), entry, i/2);
        } else {
            fmpq_zero(fmpq_mat_entry(moments, 0, i));
        }
    }

    fmpz_clear(tmp);
    fmpq_clear(entry);
    return;
}


void chebyshevu_polynomial(fmpq_poly_t Un, const int n) {
    /* Compute the n-th Chebyshev polynomial by a
     * three term recursion:
     *
     * U_0(x) = 1
     * U_1(x) = 2 * x
     *
     * U_{n+1}(x) = 2 * x * U_n(x)  -  U_{n-1}(x)
     *
     * This sequence yields Chebyshev polynomials.
     */
    fmpq_poly_t U0, U1;
    fmpq_poly_t x;
    int i;

    fmpq_poly_init(U0);
    fmpq_poly_init(U1);
    fmpq_poly_init(x);

    fmpq_poly_set_str(x, "2  0 1");
    fmpq_poly_canonicalise(x);

    fmpq_poly_one(U0);
    fmpq_poly_canonicalise(U0);

    fmpq_poly_set(U1, x);
    fmpq_poly_scalar_mul_si(U1, U1, 2);
    fmpq_poly_canonicalise(U1);

    if(n == 0) {
        fmpq_poly_set(Un, U0);
    } else if(n == 1) {
        fmpq_poly_set(Un, U1);
    } else {
        for(i = 1; i < n; i++) {
            fmpq_poly_mul(Un, x, U1);
            fmpq_poly_scalar_mul_si(Un, Un, 2);
            fmpq_poly_sub(Un, Un, U0);
            fmpq_poly_set(U0, U1);
            fmpq_poly_set(U1, Un);
        }
    }

    fmpq_poly_clear(U0);
    fmpq_poly_clear(U1);

    fmpq_poly_canonicalise(Un);
    return;
}


void integrate_chebyshevu(fmpq_t I, const int n) {
    /* Integrate
     *
     * I = \int_{-1}^{1} \sqrt{1-x^2} x^n dx
     *
     * n even:  I = \frac{\sqrt{\pi}}{2} \frac{\Gamma{\frac{n+1}{2}}}{\Gamma{2+\frac{n}{2}}}
     * n odd:   I = 0
     *
     * 1/2  0  1/8  0  1/16  0  5/128  0  7/256  0  21/1024  0  33/2048  0  429/32768
     * We omit a factor of \pi
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
        for(i = 1; i <= n/2+1; i++) {
            fmpz_set_ui(tmp, i);
            fmpq_div_fmpz(I, I, tmp);
        }
        fmpq_div_2exp(I, I, n/2 + 1);
    }

    fmpz_clear(tmp);
    return;
}


void moments_chebyshevu(fmpq_mat_t moments, const int n) {
    /*
     *
     */
    int i, j;
    fmpz_t tmp;
    fmpq_t entry;

    fmpq_mat_init(moments, 1, n);
    fmpz_init(tmp);
    fmpq_init(entry);
    fmpq_one(entry);

    for(i = 0; i < n; i++) {
        if(i % 2 == 0) {
            fmpq_one(entry);
            for(j = 1; j <= i/2; j++) {
                fmpz_set_ui(tmp, 2*j - 1);
                fmpq_mul_fmpz(entry, entry, tmp);
            }
            for(j = 1; j <= i/2+1; j++) {
                fmpz_set_ui(tmp, j);
                fmpq_div_fmpz(entry, entry, tmp);
            }
            fmpq_div_2exp(fmpq_mat_entry(moments, 0, i), entry, i/2 + 1);
        } else {
            fmpq_zero(fmpq_mat_entry(moments, 0, i));
        }
    }

    fmpz_clear(tmp);
    fmpq_clear(entry);
    return;
}


#endif
