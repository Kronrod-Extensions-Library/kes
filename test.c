/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>

#include "libkes2.h"


int main(int argc, char* argv[]) {
    int n;
    fmpq_t M;

    printf("PHY:\n");
    for(n = 0; n < 15; n++) {
	integrate_hermite_phy(M, n);
	/* M * sqrt(2pi) */
	fmpq_print(M);
	printf("  ");
    }
    printf("\n");

    printf("PRO:\n");
    for(n = 0; n < 15; n++) {
	integrate_hermite_pro(M, n);
	fmpq_print(M);
	printf("  ");
    }
    printf("\n");

    return EXIT_SUCCESS;
}
