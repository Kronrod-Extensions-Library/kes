/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>

#include "genzkeister.h"


int main(int argc, char* argv[]) {

    if(argc < 1) {
        printf("Compute Genz-Keister quadrature rule\n");
        printf("Syntax: gk [-dc D] [-dp D] [-pn] [-pw] [-pge] [-pwf] [-K K] n\n");
        printf("Options:\n");
        printf("        -dc  Compute nodes and weights up to this number of decimal digits\n");
        printf("        -dp  Print this number of decimal digits\n");
        printf("        -pn  Do not print the nodes\n");
        printf("        -pw  Do not print the weights\n");
        printf("        -pge Print the generators\n");
        printf("        -pwf Print the weight factors\n");
        printf("        -K   Set the level of the rule\n");
        return EXIT_FAILURE;
    }

    unsigned int K = 1;
    //int levels[argc-1];
    int target_prec = 53;
    int nrprintdigits = 20;
    bool print_nodes = true;
    bool print_weights = true;
    bool print_generators = false;
    bool print_weightfactors = false;

    //int k = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-dc")) {
            /* 'digits' is in base 10 and log(10)/log(2) = 3.32193 */
            target_prec = 3.32193 * atoi(argv[i+1]);
            i++;
        } else if (!strcmp(argv[i], "-dp")) {
            nrprintdigits = atoi(argv[i+1]);
            i++;
        } else if (!strcmp(argv[i], "-pn")) {
            print_nodes = false;
        } else if (!strcmp(argv[i], "-pw")) {
            print_weights = false;
        } else if (!strcmp(argv[i], "-pge")) {
            print_generators = true;
        } else if (!strcmp(argv[i], "-pwf")) {
            print_weightfactors = true;
        } else if (!strcmp(argv[i], "-K")) {
            K = atoi(argv[i+1]);
            i++;
            //} else {
            //    levels[k] = atoi(argv[i]);
            //    k++;
        }
    }

    /* Dimension of the quadrature rule */
    const unsigned int D = DIMENSION;

    // Note K is shifted by 1 wrt the original paper
    // This makes K=3 equal to the Gauss-Hermite Rule of 3 points.
    K--;

    /* Rule definition */
    std::vector<int> levels = {1, 2, 6, 10, 16};

    /* Iteratively compute quadrature nodes and weights */
    generators_t G;
    tables_t T;
    nodes_t<D> nodes;
    weights_t weights;
    rule_t<D> rule;

    for(int working_prec = target_prec; ; working_prec *= 2) {
        std::cout << "--------------------------------------------------\n";
        std::cout << "Working precision (bits): " << working_prec << std::endl;

        /* Compute data tables */
        G = compute_generators(levels, 2*working_prec);
        T = compute_weightfactors(G, 2*working_prec);

        /* Compute a Genz-Keister quadrature rule */
        rule = genz_keister_construction<D>(K, G, T, working_prec);

        nodes = rule.first;
        weights = rule.second;

        /* Accuracy goal reached? */
        if(check_accuracy<D>(nodes, weights, target_prec)) {
            break;
        }
    }

    /* Print nodes and weights */

    std::cout << "==================================================\n";
    std::cout << "DIMENSION: " << D << std::endl;
    std::cout << "LEVEL: " << K+1 << std::endl;
    std::cout << "NUMBER NODES: " << nodes.size() << std::endl;

    if(print_generators) {
        std::cout << "==================================================\n";
        std::cout << "GENERATORS" << std::endl;

        for(auto it=G.begin(); it != G.end(); it++) {
            std::cout << "| ";
            arb_printd(&(*it), nrprintdigits);
            std::cout << std::endl;
        }
    }

    if(print_weightfactors) {
        std::cout << "==================================================\n";
        std::cout << "WEIGHTFACTORS" << std::endl;

        auto WF = std::get<2>(T);
        arb_mat_printd(&WF, nrprintdigits);
    }

    if(print_nodes) {
        std::cout << "==================================================\n";
        std::cout << "NODES" << std::endl;

        for(auto it=nodes.begin(); it != nodes.end(); it++) {
            std::cout << "| ";
            for(int d=0; d < D; d++) {
                arb_printd(&((*it)[d]), nrprintdigits);
                std::cout << ",\t\t";
            }
            std::cout << std::endl;
        }
    }

    if(print_weights) {
        std::cout << "==================================================\n";
        std::cout << "WEIGHTS" << std::endl;

        for(auto it=weights.begin(); it != weights.end(); it++) {
            std::cout << "| ";
            arb_printd(&(*it), nrprintdigits);
            std::cout << std::endl;
        }
    }

    std::cout << "==================================================\n";

    return EXIT_SUCCESS;
}
