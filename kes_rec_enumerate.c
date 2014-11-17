/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "libkes.h"


int main(int argc, char* argv[]) {
    int i;
    int maxrec, n, maxp;
    fmpq_poly_t Pn;
    fmpz_mat_t table;
    int validate_weights;
    int loglevel;

    if(argc <= 1) {
        printf("Recursively search for generalized Kronrod extensions of Gauss rules\n");
        printf("Syntax: kes_rec_enumerate [-vw] [-l L] n max_p max_rec_depth\n");
        printf("        -vw  Validate the polynomial extension by weights\n");
        printf("        -l   Set the log level\n");
        return EXIT_FAILURE;
    }

    n = 1;
    maxp = 1;
    maxrec = 1;
    validate_weights = 0;
    loglevel = 0;

    for(i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-l")) {
            loglevel = atoi(argv[i+1]);
            i++;
        } else if (!strcmp(argv[i], "-vw")) {
            validate_weights = 1;
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
    printf("Search for recursive extensions of: P%i\n", n);
    printf("Maximal allowed extension order p: %i\n", maxp);
    printf("Maximal allowed recursion depth: %i\n", maxrec);
    printf("-----------------------------------------\n");

    /* Initialise the basis polynomial P1 */
    fmpq_poly_init(Pn);
    polynomial(Pn, n);

    /* Start the recursive search for E_i */
    fmpz_mat_init(table, maxrec+1, 1);
    fmpz_mat_zero(table);
    fmpz_set_ui(fmpz_mat_entry(table, 0, 0), n);

    printf("RULE: %i  ", 0);
    fmpz_print(fmpz_mat_entry(table, 0, 0));
    printf("\n");

    recursive_enumerate(Pn, maxp, 0, maxrec, table, validate_weights, loglevel);

    fmpq_poly_clear(Pn);

    return EXIT_SUCCESS;
}
