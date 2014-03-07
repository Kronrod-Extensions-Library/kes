/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "libkes2.h"


int main(int argc, char* argv[]) {
    int i, j, k;
    int levels[argc-1];
    fmpq_poly_t Pn, Ep;
    long deg;
    char *strf;
    fmpcb_ptr nodes;
    fmpcb_ptr weights;
    int solvable;
    int validate_extension, validate_weights;
    int valid;
    int comp_nodes, comp_weights;
    int target_prec;
    int nrprintdigits;
    int loglevel;

    if(argc <= 1) {
        printf("Compute a nested generalized Kronrod extension of a Gauss rule\n");
        printf("Syntax: kes [-ve] [-vw] [-cn] [-cw] [-dc D] [-dp D] [-l L] n p1 p2 ... pk\n");
        printf("Options:\n");
        printf("        -ve  Validate the polynomial extension by nodes\n");
        printf("        -vw  Validate the polynomial extension by weights\n");
        printf("        -cn  Compute the nodes\n");
        printf("        -cw  Compute the weights\n");
        printf("        -dc  Compute nodes and weights up to this number of decimal digits\n");
        printf("        -dp  Print this number of decimal digits\n");
        printf("        -l   Set the log level\n");
        return EXIT_FAILURE;
    }

    target_prec = 53;
    nrprintdigits = 20;
    loglevel = 8;
    comp_nodes = 0;
    comp_weights = 0;
    validate_extension = 0;
    validate_weights = 0;

    k = 0;
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-ve")) {
            validate_extension = 1;
        } else if (!strcmp(argv[i], "-vw")) {
            validate_weights = 1;
        } else if (!strcmp(argv[i], "-cn")) {
            comp_nodes = 1;
        } else if (!strcmp(argv[i], "-cw")) {
            comp_weights = 1;
        } else if (!strcmp(argv[i], "-dc")) {
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
            levels[k] = atoi(argv[i]);
            k++;
        }
    }

    /* Compute extension */
    fmpq_poly_init(Pn);
    polynomial(Pn, levels[0]);

    printf("Starting with polynomial:\n");
    strf = fmpq_poly_get_str_pretty(Pn, "t");
    flint_printf("P1 : %s\n", strf);
    printf("Extension levels are:");
    for(i = 0; i < k; i++) {
        printf(" %i", levels[i]);
    }
    printf("\n");

    fmpq_poly_init(Ep);
    if(k > 1) {
        solvable = find_multi_extension(Ep, Pn, k, levels, validate_extension, loglevel);
    } else {
        fmpq_poly_one(Ep);
        solvable = 1;
    }

    if(solvable) {

        fmpq_poly_mul(Pn, Pn, Ep);
        fmpq_poly_canonicalise(Pn);
        if(! fmpq_poly_is_squarefree(Pn)) {
            printf(" =====>  FINAL POLYNOMIAL NOT SQF  <=====");
        }

        if(loglevel >= 2) {
            printf("-------------------------------------------------\n");
            printf("Ending with final polynomial:\n");
            strf = fmpq_poly_get_str_pretty(Pn, "t");
            flint_printf("P : %s\n", strf);
        }

        deg = fmpq_poly_degree(Pn);
        nodes = _fmpcb_vec_init(deg);
        weights = _fmpcb_vec_init(deg);

        if(comp_weights || validate_weights) {
            compute_nodes_and_weights(nodes, weights, Pn, target_prec, loglevel);
        } else if(comp_nodes) {
            compute_nodes(nodes, Pn, target_prec, loglevel);
        }

        valid = 1;
        if(validate_weights) {
            valid = validate_extension_by_weights(weights, deg, target_prec, loglevel);
        }

        if(! valid) {
            printf("**************************************\n");
            printf("*** EXTENSION WITH INVALID WEIGHTS ***\n");
            printf("**************************************\n");
        }

        /* Print roots and weights */
        if((validate_weights && valid) || ! validate_weights) {
            if(comp_nodes) {
                printf("-------------------------------------------------\n");
                printf("The nodes are:\n");
                for(j = 0; j < deg; j++) {
                    printf("| ");
                    fmpcb_printd(nodes + j, nrprintdigits);
                    printf("\n");
                }
            }

            if(comp_weights) {
                printf("-------------------------------------------------\n");
                printf("The weights are:\n");
                for(j = 0; j < deg; j++) {
                    printf("| ");
                    fmpcb_printd(weights + j, nrprintdigits);
                    printf("\n");
                }
            }
        }

        _fmpcb_vec_clear(nodes, deg);
        _fmpcb_vec_clear(weights, deg);
    }

    flint_free(strf);
    fmpq_poly_clear(Pn);
    fmpq_poly_clear(Ep);
    return EXIT_SUCCESS;
}
