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

    if(argc < 2) {
        printf("Compute Genz-Keister quadrature rule\n");
        printf("Syntax: genzkeister [-dc D] [-dp D] [-pn] [-pw] [-pge] [-pwf] -K K [n1 n2 n3 ...nk]\n");
        printf("Options:\n");
        printf("        -dc  Compute nodes and weights up to this number of decimal digits\n");
        printf("        -dp  Print this number of decimal digits\n");
        printf("        -pn  Do not print the nodes\n");
        printf("        -pw  Do not print the weights\n");
        printf("        -pge Print the generators\n");
        printf("        -pmo Print the moments\n");
        printf("        -pai Print the a_i values\n");
        printf("        -pwf Print the weight factors\n");
        printf("        -pzs Print the Z-sequence\n");
        printf("        -K   Set the level of the rule\n");
        printf("        Further optional parameters  n1 ... nk  define the Kronrod extension used\n");
        return EXIT_FAILURE;
    }

    unsigned int K = 1;
    int target_prec = 53;
    int nrprintdigits = 20;
    std::vector<int> levels;
    bool print_nodes = true;
    bool print_weights = true;
    bool print_generators = false;
    bool print_moments = false;
    bool print_avalues = false;
    bool print_weightfactors = false;
    bool print_zsequence = false;

    for(int i = 1; i < argc; i++) {
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
        } else if (!strcmp(argv[i], "-pmo")) {
            print_moments = true;
        } else if (!strcmp(argv[i], "-pai")) {
            print_avalues = true;
        } else if (!strcmp(argv[i], "-pwf")) {
            print_weightfactors = true;
        } else if (!strcmp(argv[i], "-pzs")) {
            print_zsequence = true;
        } else if (!strcmp(argv[i], "-K")) {
            K = atoi(argv[i+1]);
            i++;
        } else {
            levels.push_back(atoi(argv[i]));
        }
    }

    /* Dimension of the quadrature rule */
    const unsigned int D = DIMENSION;

    // Note K is shifted by 1 wrt the original paper
    // This makes K=3 equal to the Gauss Rule of 3 points.
    K--;

    /* Default rule definition */
    if(levels.size() == 0) {
#ifdef LEGENDRE
        levels = {1, 2, 4, 8, 16, 32};//, 64};
#endif
#ifdef HERMITEPRO
        levels = {1, 2, 6, 10, 16};//, 68};
#endif
#ifdef HERMITE
        levels = {1, 2, 6, 10, 16};//, 68};
#endif
#ifdef CHEBYSHEVT
        levels = {1, 2, 4, 6, 12, 24};//, 48};
#endif
#ifdef CHEBYSHEVU
        levels = {1, 2, 4, 8, 16, 32};//, 64};
#endif
    }

    std::cout << "Kronrod extension:  ";
    for(auto it=levels.begin(); it != levels.end(); it++) {
        std::cout << *it << "  ";
    }
    std::cout << "\n==================================================" << std::endl;

    /* Iteratively compute quadrature nodes and weights */
    generators_t G;
    tables_t T;
    nodes_t<D> nodes;
    weights_t weights;
    rule_t<D> rule;

    for(unsigned int working_prec = target_prec; ; working_prec *= 2) {
        std::cout << "--------------------------------------------------\n";
        std::cout << "Working precision (bits): " << working_prec << std::endl;

        /* Compute data tables */
        // TODO: Assert generator g_0 = 0
        G = compute_generators(levels, 2*working_prec);
        T = compute_tables(G, 2*working_prec);

        if(K >= G.size()) {
            std::cout << "***********************************\n";
            std::cout << "*** not enough generators found ***\n";
            std::cout << "***********************************\n";
            break;
        }

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

    if(print_moments) {
        std::cout << "==================================================\n";
        std::cout << "MOMENTS" << std::endl;

        moments_t M = std::get<0>(T);
        fmpq_mat_print(&M);
        std::cout << std::endl;
    }

    if(print_avalues) {
        std::cout << "==================================================\n";
        std::cout << "A-VALUES" << std::endl;

        ai_t A = std::get<1>(T);
        arb_mat_printd(&A, nrprintdigits);
    }

    if(print_weightfactors) {
        std::cout << "==================================================\n";
        std::cout << "WEIGHTFACTORS" << std::endl;

        wft_t WF = std::get<2>(T);
        arb_mat_printd(&WF, nrprintdigits);
    }

    if(print_zsequence) {
        std::cout << "==================================================\n";
        std::cout << "Z-SEQUENCE" << std::endl;

        z_t Z = std::get<3>(T);
        std::cout << "Z = [ ";
        for(auto it=Z.begin(); it != Z.end(); it++) {
            std::cout << (*it) << " ";
        }
        std::cout << "]" << std::endl;
    }

    if(print_nodes) {
        std::cout << "==================================================\n";
        std::cout << "NODES" << std::endl;

        for(auto it=nodes.begin(); it != nodes.end(); it++) {
            std::cout << "| ";
            for(unsigned int d=0; d < D; d++) {
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
