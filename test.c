/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>

#include "polynomials.h"
#include "libkes2.h"


int main(int argc, char* argv[]) {
    int n;
    fmpq_t M;
    fmpq_poly_t P;
    fmpq_poly_t L;
    fmpq_poly_t H;
    char *strf;

    fmpq_poly_init(P);
    fmpq_poly_init(L);
    fmpq_poly_init(H);

    printf("Legendre Polynomials:\n");
    for(n = 0; n < 10; n++) {
	legendre_polynomial(P, n);
	strf = fmpq_poly_get_str_pretty(P, "t");
	flint_printf("L%d : %s\n", n, strf);
    }
    printf("\n");

    printf("Legendre Moments:\n");
    for(n = 0; n < 15; n++) {
    	integrate_legendre(M, n);
    	fmpq_print(M);
    	printf("  ");
    }
    printf("\n");


    printf("Laguerre Polynomials:\n");
    for(n = 0; n < 10; n++) {
    	laguerre_polynomial(L, n);
    	strf = fmpq_poly_get_str_pretty(L, "t");
    	flint_printf("L%d : %s\n", n, strf);
    }
    printf("\n");

    printf("Laguerre Moments:\n");
    for(n = 0; n < 15; n++) {
    	integrate_laguerre(M, n);
    	fmpq_print(M);
    	printf("  ");
    }
    printf("\n");


    printf("Hermite PHY:\n");
    for(n = 0; n < 10; n++) {
    	hermite_polynomial_phy(H, n);
    	strf = fmpq_poly_get_str_pretty(H, "t");
    	flint_printf("H%d : %s\n", n, strf);
    }
    printf("\n");

    printf("Hermite PHY Moments:\n");
    for(n = 0; n < 15; n++) {
    	integrate_hermite_phy(M, n);
    	fmpq_print(M);
    	printf("  ");
    }
    printf("\n");


    printf("Hermite PRO:\n");
    for(n = 0; n < 10; n++) {
    	hermite_polynomial_pro(H, n);
    	strf = fmpq_poly_get_str_pretty(H, "t");
    	flint_printf("H%d : %s\n", n, strf);
    }
    printf("\n");

    printf("Hermite PRO Moments:\n");
    for(n = 0; n < 15; n++) {
    	integrate_hermite_pro(M, n);
    	fmpq_print(M);
    	printf("  ");
    }
    printf("\n");


    return EXIT_SUCCESS;
}
