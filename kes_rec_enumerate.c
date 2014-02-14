/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "libkes2.h"


int main(int argc, char* argv[]) {
    int i;
    int maxrec, n, maxp;
    fmpq_poly_t Hn;
    fmpz_mat_t table;
    int loglevel;

    if(argc <= 1) {
        printf("Recursively search for generalized Kronrod extensions of Gauss-Hermite rules\n");
        printf("Syntax: kes_rec_enumerate [-l L] n max_p max_rec_depth\n");
        printf("        -l   Set the log level\n");
        return EXIT_FAILURE;
    }

    n = 1;
    maxp = 1;
    maxrec = 1;
    loglevel = 0;

    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-l")) {
            loglevel = atoi(argv[i+1]);
            i++;
        } else {
            if(i + 2 < argc) {
                n = atoi(argv[i]);
                maxp = atoi(argv[i+1]);
                maxrec = atoi(argv[i+2]);
                break;
            } else {
                printf("Missing values for n or max_p!\n");
                return EXIT_FAILURE;
            }
        }
    }

    printf("-----------------------------------------\n");
    printf("Search for recursive extensions of: H%i\n", n);
    printf("Maximal allowed extension order p: %i\n", maxp);
    printf("Maximal allowed recursion depth: %i\n", maxrec);
    printf("-----------------------------------------\n");

    /* Initialise the basis polynomial P1 */
    fmpq_poly_init(Hn);
    hermite_polynomial_pro(Hn, n);

    /* Start the recursive search for E_i */
    fmpz_mat_init(table, maxrec+3, 1);
    fmpz_mat_zero(table);
    fmpz_set_ui(fmpz_mat_entry(table, 0, 0), n);

    printf("RULE: %i  ", 1);
    fmpz_mat_print(table);
    printf("\n");

    recursive_enumerate(Hn, n, maxp, 0, maxrec, table, loglevel);

    fmpq_poly_clear(Hn);

    return EXIT_SUCCESS;
}
