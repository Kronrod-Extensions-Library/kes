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

//#include "boost/multi_array.hpp"

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
    long deg;
    acb_ptr generators;

    int target_prec = 530;
    //int zero_eps = 53;
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

    int NUMGEN = G.size();


    /* Compute moments of exp(-x^2/2) */
    fmpz_mat_t M;
    fmpz_mat_init(M, 1, 2*NUMGEN+1);
    fmpz_t I, tmp;
    fmpz_init(I);
    fmpz_one(I);
    for(long i=0; i < 2*NUMGEN+1; i++) {
	if(i % 2 == 0) {
	    fmpz_set(fmpz_mat_entry(M, 0, i), I);
            fmpz_set_ui(tmp, i + 1);
            fmpz_mul(I, I, tmp);
	} else {
	    fmpz_set_ui(fmpz_mat_entry(M, 0, i), 0);
	}
    }

    std::cout << "=========================================\n";
    fmpz_mat_print_pretty(M);


    /* Compute the values of all a_i */
    arb_mat_t A;
    arb_mat_init(A, 1, NUMGEN+1);

    arb_poly_t term, poly;
    arb_poly_init(term);
    arb_poly_init(poly);
    arb_poly_one(poly);

    arb_t t, u, c, ai;
    arb_init(t);
    arb_init(u);
    arb_init(c);
    arb_init(ai);

    // a_0 = 1
    arb_set_ui(arb_mat_entry(A, 0, 0), 1);

    // a_i
    int i = 1;
    for(auto it=G.begin(); it != G.end(); it++) {
	// Construct the polynomial term by term
	arb_poly_set_coeff_si(term, 2, 1);
	arb_pow_ui(t, &(*it), 2, target_prec);
	arb_neg(t, t);
	arb_poly_set_coeff_arb(term, 0, t);
	arb_poly_mul(poly, poly, term, target_prec);
	// Compute the contributions to a_i for each monomial
	arb_zero(ai);
	long deg = arb_poly_degree(poly);
	for(long d=0; d <= deg; d++) {
	    arb_set_fmpz(t, fmpz_mat_entry(M, 0, d));
	    arb_poly_get_coeff_arb(c, poly, d);
	    arb_mul(t, c, t, target_prec);
	    arb_add(ai, ai, t, target_prec);
	}
	// TODO: More zero tests
	if(arf_is_zero(arb_midref(ai))) {
	    arb_zero(ai);
	}
	// Assign a_i
	arb_set(arb_mat_entry(A, 0, i), ai);
	i++;
    }

    std::cout << "=========================================\n";
    arb_mat_printd(A, 5);


    /* Precompute all weight factors */
    arb_mat_t WF;
    arb_mat_init(WF, NUMGEN+1, NUMGEN+1);
    arb_mat_zero(WF);

    for(int xi=0; xi <= NUMGEN; xi++) {
	arb_one(c);
	for(int theta=0; theta <= NUMGEN; theta++) {
	    if(theta != xi) {
		arb_pow_ui(t, &G[theta], 2, target_prec);
		arb_pow_ui(u, &G[xi], 2, target_prec);
		arb_sub(t, u, t, target_prec);
		arb_mul(c, c, t, target_prec);
	    }
	    if(theta >= xi) {
		arb_div(t, arb_mat_entry(A, 0, theta), c, target_prec);
		arb_set(arb_mat_entry(WF, xi, theta), t);
	    }
	}
    }

    std::cout << "=========================================\n";
    arb_mat_printd(WF, 5);







    fmpq_poly_clear(Pn);
    fmpq_poly_clear(Ep);

    return EXIT_SUCCESS;
}
