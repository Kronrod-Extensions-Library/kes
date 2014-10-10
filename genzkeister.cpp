/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <list>
#include <vector>
#include <iostream>

#include "arf.h"
#include "arb.h"
#include "acb.h"

#include "libkes2.h"


void maxminsort(std::vector<arb_struct>& G, acb_ptr g, long n) {
    // Put all elements in t
    std::list<arb_struct> t(0);
    for(int i = 0; i < n; i++) {
	if(arb_is_nonnegative(acb_realref(g+i))) {
	    t.push_back(acb_realref(g+i)[0]);
	}
    }
    // Sort them by radius, drop negative ones
    bool largest = true;
    while(t.size() > 0) {
      	auto I = t.begin();
	// Look for largest or smallest element
	for(auto it=t.begin(); it != t.end(); it++) {
	    if(largest) {
		if(arf_cmp(arb_midref(&(*I)), arb_midref(&(*it))) <= 0) {
		    I = it;
		}
	    } else {
		if(arf_cmp(arb_midref(&(*I)), arb_midref(&(*it))) >= 0) {
		    I = it;
		}
	    }
	}
	// Copy over into G
	G.push_back(*I);
	t.erase(I);
	largest = !largest;
    }
}



int main(int argc, char* argv[]) {
    fmpq_poly_t Pn;
    fmpq_poly_t Ep;
    char *strf;
    long deg;
    acb_ptr generators;

    int target_prec = 53;
    int nrprintdigits = 20;
    int loglevel = 8;
    int validate_extension = 0;


    if(argc <= 1) {
        printf("Compute Genz-Keister quadrature rule\n");
        printf("Syntax: gk [-dc D] [-dp D] [-l L] n\n");
        printf("Options:\n");
        printf("        -dc  Compute nodes and weights up to this number of decimal digits\n");
        printf("        -dp  Print this number of decimal digits\n");
        printf("        -l   Set the log level\n");
        return EXIT_FAILURE;
    }

    /* Rule definition */
    const int k = 5;
    int levels[k] = {1, 2, 6, 10, 16};

    /* Compute generators */
    std::vector<arb_struct> G(0);

    /* Compute extension recursively */
    fmpq_poly_init(Pn);
    fmpq_poly_init(Ep);

    polynomial(Pn, levels[0]);

    deg = fmpq_poly_degree(Pn);
    generators = _acb_vec_init(deg);
    compute_nodes(generators, Pn, target_prec, loglevel);
    maxminsort(G, generators, deg);
    //_acb_vec_clear(generators, deg);

    for(int i = 1; i < k; i++) {
	bool solvable = find_extension(Ep, Pn, levels[i], loglevel);
	if(!solvable) {
	    break;
	}

	/* Compute generators (zeros of Pn) */
	deg = fmpq_poly_degree(Ep);
	generators = _acb_vec_init(deg);
	compute_nodes(generators, Ep, target_prec, loglevel);
	maxminsort(G, generators, deg);
	//_acb_vec_clear(generators, deg);

        /* Iterate */
        fmpq_poly_mul(Pn, Pn, Ep);
        fmpq_poly_canonicalise(Pn);
    }

    // Print
    std::cout << "=========================================\n";
    for(auto it=G.begin(); it != G.end(); it++) {
	std::cout << "||  ";
        arb_printd(&(*it), nrprintdigits);
        std::cout << std::endl;
    }

    flint_free(strf);
    fmpq_poly_clear(Pn);
    fmpq_poly_clear(Ep);

    return EXIT_SUCCESS;
}
