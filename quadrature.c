/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "libkes2.h"


int main(int argc, char* argv[]) {
    int i, j;
    int deg;
    fmpq_poly_t Pn;
    char *strf;
    acb_ptr nodes;
    acb_ptr weights;
    int working_prec;
    int target_prec;
    int nrprintdigits;
    int loglevel;

    if(argc <= 1) {
        printf("Compute Gauss quadrature rule\n");
        printf("Syntax: quadrature [-dc D] [-dp D] [-l L] n\n");
        printf("Options:\n");
        printf("        -dc  Compute nodes and weights up to this number of decimal digits\n");
        printf("        -dp  Print this number of decimal digits\n");
        printf("        -l   Set the log level\n");
        return EXIT_FAILURE;
    }

    deg = 0;
    target_prec = 53;
    nrprintdigits = 20;
    loglevel = 8;

    for(i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-dc")) {
            /* 'digits' is in base 10 and log(10)/log(2) = 3.32193 */
            target_prec = 3.32193 * atoi(argv[i+1]);
            i++;
        } else if (!strcmp(argv[i], "-dp")) {
            nrprintdigits = atoi(argv[i+1]);
            i++;
        } else if (!strcmp(argv[i], "-l")) {
            loglevel = atoi(argv[i+1]);
            i++;
        } else {
            deg = atoi(argv[i]);
        }
    }

    /* Compute polynomial */
    fmpq_poly_init(Pn);
    polynomial(Pn, deg);

    printf("Starting with polynomial:\n");
    strf = fmpq_poly_get_str_pretty(Pn, "t");
    flint_printf("P : %s\n", strf);

    nodes = _acb_vec_init(deg);
    weights = _acb_vec_init(deg);

    /* Precision in number of bits */
    for(working_prec = target_prec; ; working_prec *= 2) {
        /* Find nodes and weights */
        compute_nodes(nodes, Pn, working_prec, loglevel);
        sort_nodes(nodes, deg);
        evaluate_weights_formula(weights, nodes, deg, working_prec);
        /* Accuracy goal reached? */
        if(check_accuracy(nodes, deg, target_prec) && check_accuracy(weights, deg, target_prec)) {
            break;
        }
    }

    /* Print roots and weights */
    printf("-------------------------------------------------\n");
    printf("The nodes are:\n");
    for(j = 0; j < deg; j++) {
        printf("| ");
        acb_printd(nodes + j, nrprintdigits);
        printf("\n");
    }
    printf("-------------------------------------------------\n");
    printf("The weights are:\n");
    for(j = 0; j < deg; j++) {
        printf("| ");
        acb_printd(weights + j, nrprintdigits);
        printf("\n");
    }

    _acb_vec_clear(nodes, deg);
    _acb_vec_clear(weights, deg);

    flint_free(strf);
    fmpq_poly_clear(Pn);
    return EXIT_SUCCESS;
}
