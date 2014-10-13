/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <list>
#include <vector>
#include <array>
#include <iostream>
#include <utility>

#include "arf.h"
#include "arb.h"
#include "acb.h"

#include "libkes2.h"
#include "enumerators.h"


void maxminsort(std::vector<arb_struct>& generators, acb_ptr g, long n) {
    /* Sort generators according to the max-min heuristic.
     */
    std::vector<arb_struct> t(0);
    for(int i = 0; i < n; i++) {
        if(arb_is_nonnegative(acb_realref(g+i))) {
            t.push_back(acb_realref(g+i)[0]);
        }
    }
    // Sort balls by midpoint, drop negative ones
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
        // Copy over into generators
        generators.push_back(*I);
        t.erase(I);
        largest = !largest;
    }
}


std::vector<arb_struct> compute_generators(const int k,
                                           const int levels[],
                                           const int working_prec) {
    /* Compute the generators.
     *
     * k: The number of nested extensions in the tower
     * levels: An array with the extension levels p_0, ..., p_{k-1}
     * working_prec: Working precision
     */
    std::vector<arb_struct> G(0);

    /* Compute extension recursively and obtain generators */
    fmpq_poly_t Pn;
    fmpq_poly_t Ep;
    fmpq_poly_init(Pn);
    fmpq_poly_init(Ep);

    hermite_polynomial_pro(Pn, levels[0]);

    long deg = fmpq_poly_degree(Pn);
    acb_ptr generators = _acb_vec_init(deg);
    compute_nodes(generators, Pn, working_prec, 0);
    maxminsort(G, generators, deg);
    //_acb_vec_clear(generators, deg);

    for(int i = 1; i < k; i++) {
        bool solvable = find_extension(Ep, Pn, levels[i], 0);
        if(!solvable) {
            break;
        }

        /* Compute generators (zeros of Pn) */
        deg = fmpq_poly_degree(Ep);
        generators = _acb_vec_init(deg);
        compute_nodes(generators, Ep, working_prec, 0);
        maxminsort(G, generators, deg);
        //_acb_vec_clear(generators, deg);

        /* Iterate */
        fmpq_poly_mul(Pn, Pn, Ep);
        fmpq_poly_canonicalise(Pn);
    }

    fmpq_poly_clear(Pn);
    fmpq_poly_clear(Ep);

    return G;
}


arb_mat_struct compute_weightfactors(const std::vector<arb_struct>& generators,
                                     const int working_prec) {
    /* Compute the weight factors.
     *
     * generators: List of generators
     * working_prec: Working precision
     */
    int number_generators = generators.size();

    /* Compute moments of exp(-x^2/2) */
    fmpz_mat_t M;
    fmpz_mat_init(M, 1, 2*number_generators+1);
    fmpz_t I, tmp;
    fmpz_init(I);
    fmpz_one(I);
    for(long i=0; i < 2*number_generators+1; i++) {
        if(i % 2 == 0) {
            fmpz_set(fmpz_mat_entry(M, 0, i), I);
            fmpz_set_ui(tmp, i + 1);
            fmpz_mul(I, I, tmp);
        } else {
            fmpz_set_ui(fmpz_mat_entry(M, 0, i), 0);
        }
    }

    /* Compute the values of all a_i */
    arb_mat_t A;
    arb_mat_init(A, 1, number_generators+1);

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
    for(auto it=generators.begin(); it != generators.end(); it++) {
        // Construct the polynomial term by term
        arb_poly_set_coeff_si(term, 2, 1);
        arb_pow_ui(t, &(*it), 2, working_prec);
        arb_neg(t, t);
        arb_poly_set_coeff_arb(term, 0, t);
        arb_poly_mul(poly, poly, term, working_prec);
        // Compute the contributions to a_i for each monomial
        arb_zero(ai);
        long deg = arb_poly_degree(poly);
        for(long d=0; d <= deg; d++) {
            arb_set_fmpz(t, fmpz_mat_entry(M, 0, d));
            arb_poly_get_coeff_arb(c, poly, d);
            arb_mul(t, c, t, working_prec);
            arb_add(ai, ai, t, working_prec);
        }
        // TODO: More zero tests
        if(arf_is_zero(arb_midref(ai))) {
            arb_zero(ai);
        }
        // Assign a_i
        arb_set(arb_mat_entry(A, 0, i), ai);
        i++;
    }

    /* Precompute all weight factors */
    arb_mat_t weight_factors;
    arb_mat_init(weight_factors, number_generators+1, number_generators+1);
    arb_mat_zero(weight_factors);

    for(int xi=0; xi <= number_generators; xi++) {
        arb_one(c);
        for(int theta=0; theta <= number_generators; theta++) {
            if(theta != xi) {
                arb_pow_ui(t, &generators[theta], 2, working_prec);
                arb_pow_ui(u, &generators[xi], 2, working_prec);
                arb_sub(t, u, t, working_prec);
                arb_mul(c, c, t, working_prec);
            }
            if(theta >= xi) {
                arb_div(t, arb_mat_entry(A, 0, theta), c, working_prec);
                arb_set(arb_mat_entry(weight_factors, xi, theta), t);
            }
        }
    }

    return *weight_factors;
}


template <int D>
std::vector<std::array<arb_struct, D>>
Points(const std::array<int, D> P,
       const std::vector<arb_struct>& generators,
       const int working_prec) {
    /* Compute fully symmetric quadrature nodes for given partition `P`.
     *
     * P: Partition `P`
     * generators: Table with precomputed generators
     * working_prec: Working precision
     */
    std::vector<std::array<arb_struct, D>> points;
    std::array<int, D> Q;

    // Number of 0 entries in P
    int xi = nz<D>(P);
    int nsf = xi<D ? (2<<(D-xi-1)) : 1;

    std::list<std::array<int, D>> permutations = Permutations<D>(P);
    for(auto it=permutations.begin(); it != permutations.end(); it++) {
        Q = *it;
        for(int v=0; v < nsf; v++) {
            std::array<arb_struct, D> point;
            int u = 0;
            for(int d=0; d < D; d++) {
                arb_t t;
                arb_init(t);
                // Generators corresponding to partition Q give current node
                int i = Q[d];
                arb_set(t, &(generators[i]));
                // Compute sign flip
                if(i != 0) {
                    if((v >> u) & 1) {
                        arb_neg(t, t);
                    }
                    u++;
                }
                // Assign
                point[d] = *t;
            }
            points.push_back(point);
        }
    }
    return points;
}


template <int D>
std::vector<arb_struct>
Weights(std::array<int, D> P,
        int K,
        arb_mat_struct& weight_factors,
        int working_prec) {
    /* Function to compute weights for partition `P`.
     *
     * P: Partition `P`
     * K: Rule order
     * weight_factors: Table with precomputed weight factors
     * working_prec: Working precision
     */
    std::vector<arb_struct> weights;

    arb_t W, w;
    arb_init(W);
    arb_init(w);
    arb_zero(W);

    std::list<std::array<int, D>> U = LatticePoints<D>(K-sum<D>(P));
    std::array<int, D> Q;

    for(auto it=U.begin(); it != U.end(); it++) {
        Q = *it;
        arb_one(w);
        for(int d=0; d < D; d++) {
            arb_mul(w, w, arb_mat_entry(&weight_factors, P[d], P[d]+Q[d]), working_prec);
        }
        arb_add(W, W, w, working_prec);
    }
    // Number of non-zero
    int k = nnz<D>(P);
    if(k > 0) {
        arb_div_ui(W, W, (2 << (k-1)) , working_prec);
    }
    weights.push_back(*W);

    arb_clear(w);
    //arb_clear(W);

    return weights;
}


template <int D>
std::pair<std::vector<std::array<arb_struct, D> >, std::vector<arb_struct> >
GenzKeisterConstruction(int K,
                        const std::vector<arb_struct>& generators,
                        arb_mat_struct& weight_factors,
                        const int Z[],
                        int working_prec) {
    /* Compute the Genz-Keister construction
     *
     * K: Level of the quadrature rule
     * generators: Table with precomputed generators
     * weight_factors: Table with the weight factors
     * Z:
     * working_prec: Working precision
     */
    std::vector<std::array<arb_struct, D>> points;
    std::vector<arb_struct> weights;

    // Iterate over all relevant integer partitions
    std::list<std::array<int, D>> partitions = Partitions<D>(K);
    for(auto it=partitions.begin(); it != partitions.end(); it++) {
        //
        std::array<int, D> P = *it;
        int s = 0;
        for(int d=0; d < D; d++) {
            s += P[d];
            s += Z[P[d]];
        }
        //
        if(s <= K) {
            // Compute nodes and weights for given partition
            std::vector<std::array<arb_struct, D>> p = Points<D>(P, generators, working_prec);
            std::vector<arb_struct> w = Weights<D>(P, K, weight_factors, working_prec);
            for(int i=0; i < p.size(); i++) {
                points.push_back(p[i]);
                weights.push_back(w[0]);
            }
        }
    }

    return std::make_pair(points, weights);
}


int main(int argc, char* argv[]) {
    // Some constants
    int target_prec = 53;
    int working_prec = target_prec;

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

    const int D = 2;
    int K = 3;
    K -= 1;

    /* Rule definition */

    const int k = 5;
    int levels[k] = {1, 2, 6, 10, 16};

    /* Precompute the Z-sequence */
    // TODO Based on formula
    int Z[27] = {0,0,
                 1,0,0,
                 3,2,1,0,0,
                 5,4,3,2,1,0,0,0,
                 8,7,6,5,4,3,2,1,0
    };

    /* Compute data tables */

    std::vector<arb_struct> G = compute_generators(k, levels, working_prec);
    arb_mat_struct WF = compute_weightfactors(G, working_prec);

    /* Compute a Genz-Keister quadrature rule */

    std::pair<std::vector<std::array<arb_struct, D>>, std::vector<arb_struct>> rule;
    rule = GenzKeisterConstruction<D>(K, G, WF, Z, working_prec);

    std::vector<std::array<arb_struct, D>> points;
    std::vector<arb_struct> weights;
    points = rule.first;
    weights = rule.second;

    /* Print nodes and weights */

    std::cout << "==================================================\n";
    std::cout << "GENERATORS" << std::endl;

    for(auto it=G.begin(); it != G.end(); it++) {
        std::cout << "| ";
        arb_printd(&(*it), nrprintdigits);
        std::cout << std::endl;
    }

    std::cout << "==================================================\n";
    std::cout << "WEIGHTFACTORS" << std::endl;

    arb_mat_printd(&WF, 5);

    std::cout << "==================================================\n";
    std::cout << "NODES" << std::endl;

    for(auto it=points.begin(); it != points.end(); it++) {
        std::cout << "| ";
        for(int d=0; d < D; d++) {
            arb_printd(&((*it)[d]), 7);
            std::cout << ",\t\t";
        }
        std::cout << std::endl;
    }

    std::cout << "==================================================\n";
    std::cout << "WEIGHTS" << std::endl;

    for(auto it=weights.begin(); it != weights.end(); it++) {
        std::cout << "| ";
        arb_printd(&(*it), nrprintdigits);
        std::cout << std::endl;
    }

    std::cout << "==================================================\n";

    return EXIT_SUCCESS;
}
