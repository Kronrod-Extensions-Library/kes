/*  Author: R. Bourquin
 *  Copyright: (C) 2014 R. Bourquin
 *  License: GNU GPL v2 or above
 *
 *  A library of helper functions to search for
 *  Kronrod extensions of Gauss quadrature rules.
 */

#ifndef __HH__genzkeister
#define __HH__genzkeister

#include <list>
#include <vector>
#include <array>
#include <utility>
#include <tuple>
#include <iostream>

#include "arf.h"
#include "arb.h"
#include "acb.h"

#include "libkes2.h"
#include "enumerators.h"


typedef std::vector<arb_struct> generators_t;
typedef fmpq_mat_struct moments_t;
typedef arb_mat_struct ai_t;
typedef std::vector<int> z_t;
typedef arb_mat_struct wft_t;
typedef std::tuple<moments_t, ai_t, wft_t, z_t> tables_t;
typedef std::vector<arb_struct> weights_t;
template<int D> using node_t = std::array<arb_struct, D>;
template<int D> using nodes_t = std::vector<node_t<D>>;
template<int D> using rule_t = std::pair<nodes_t<D>, weights_t>;


void maxminsort(generators_t& generators, const acb_ptr g, const long n) {
    /* Sort generators according to the max-min heuristic.
     */
    generators_t t(0);
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


generators_t compute_generators(const std::vector<int> levels,
                                const int working_prec) {
    /* Compute the generators.
     *
     * levels: An array with the extension levels p_0, ..., p_{k-1}
     * working_prec: Working precision
     */
    int prec = 2 * working_prec;
    generators_t G(0);

    /* Compute extension recursively and obtain generators */
    fmpq_poly_t Pn;
    fmpq_poly_t Ep;
    fmpq_poly_init(Pn);
    fmpq_poly_init(Ep);

    polynomial(Pn, levels[0]);

    long deg = fmpq_poly_degree(Pn);
    acb_ptr generators = _acb_vec_init(deg);
    compute_nodes(generators, Pn, prec, 0);
    maxminsort(G, generators, deg);
    //_acb_vec_clear(generators, deg);

    for(unsigned int i = 1; i < levels.size(); i++) {
        bool solvable = find_extension(Ep, Pn, levels[i], 0);
        if(!solvable) {
            std::cout << "******************************\n";
            std::cout << "*** EXTENSION NOT SOVLABLE ***\n";
            std::cout << "******************************\n";
            break;
        }

        /* Compute generators (zeros of Pn) */
        deg = fmpq_poly_degree(Ep);
        generators = _acb_vec_init(deg);
        compute_nodes(generators, Ep, prec, 0);
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


moments_t
compute_moments(const int N) {
    /* Compute the moments of x^n with exp(-x^2/2)
     *
     * N: The number of moments to be computed
     */
    fmpq_mat_t M;
    moments(M, N);
    return *M;
}


ai_t
compute_ai(const generators_t& generators,
           const moments_t& moments,
           const int working_prec) {
    /* Compute the values of all a_i
     *
     * generators: Table with precomputed generators
     * moments: Table with precomputed moments of x^n
     * working_prec: Working precision
     */
    int number_generators = generators.size();

    arb_mat_t A;
    arb_mat_init(A, 1, number_generators+1);

    arb_poly_t term, poly;
    arb_poly_init(term);
    arb_poly_init(poly);
    arb_poly_one(poly);

    arb_t c, t, ai;
    arb_init(c);
    arb_init(t);
    arb_init(ai);

    // a_0
    arb_set_fmpq(arb_mat_entry(A, 0, 0), fmpq_mat_entry(&moments, 0, 0), working_prec);

    // a_i
    int i = 1;
    for(auto it=generators.begin(); it != generators.end(); it++) {
        // Construct the polynomial term by term
        arb_pow_ui(t, &(*it), 2, working_prec);
        arb_neg(t, t);
        arb_poly_set_coeff_si(term, 2, 1);
        arb_poly_set_coeff_arb(term, 0, t);
        arb_poly_mul(poly, poly, term, working_prec);
        // Compute the contributions to a_i for each monomial
        arb_zero(ai);
        long deg = arb_poly_degree(poly);
        for(long d=0; d <= deg; d++) {
            arb_set_fmpq(t, fmpq_mat_entry(&moments, 0, d), working_prec);
            arb_poly_get_coeff_arb(c, poly, d);
            arb_mul(t, c, t, working_prec);
            arb_add(ai, ai, t, working_prec);
        }
        // Zero tests for a_i
        if(arf_cmpabs_2exp_si(arb_midref(ai), -working_prec/2) < 0) {
            arb_zero(ai);
        }
        // Assign a_i
        arb_set(arb_mat_entry(A, 0, i), ai);
        i++;
    }

    return *A;
}


wft_t
compute_weightfactors(const generators_t& generators,
                      const ai_t& A,
                      const int working_prec) {
    /* Precompute all weight factors
     *
     * generators: Table with precomputed generators
     * A: Table with precomputed a_i values
     * working_prec: Working precision
     */
    int number_generators = generators.size();

    arb_t c, t, u;
    arb_init(t);
    arb_init(u);
    arb_init(c);

    arb_mat_t weight_factors;
    arb_mat_init(weight_factors, number_generators, number_generators);
    arb_mat_zero(weight_factors);

    for(int xi=0; xi < number_generators; xi++) {
        arb_one(c);
        for(int theta=0; theta < number_generators; theta++) {
            if(theta != xi) {
                arb_pow_ui(t, &generators[theta], 2, working_prec);
                arb_pow_ui(u, &generators[xi], 2, working_prec);
                arb_sub(t, u, t, working_prec);
                arb_mul(c, c, t, working_prec);
            }
            if(theta >= xi) {
                arb_div(t, arb_mat_entry(&A, 0, theta), c, working_prec);
                arb_set(arb_mat_entry(weight_factors, xi, theta), t);
            }
        }
    }

    return *weight_factors;
}


z_t
compute_z_sequence(const ai_t& A) {
    /* Precompute the Z-sequence
     *
     * A: Table with precomputed a_i values
     */
    int n = arb_mat_ncols(&A);

    std::vector<int> Z(0);

    int v = 0;
    for(int i=0; i < n; i++) {
        if(v == 0) {
            while( arb_is_zero(arb_mat_entry(&A, 0, i+v)) ) {
                v++;
            }
        } else {
            v--;
        }
        Z.push_back(v);
    }

    /* // TODO Based on formula */
    /* std::vector<int> Z = {0,0, */
    /*                       1,0,0, */
    /*                       3,2,1,0,0, */
    /*                       5,4,3,2,1,0,0,0, */
    /*                       8,7,6,5,4,3,2,1,0}; */

    return Z;
}


tables_t
compute_tables(const generators_t& generators,
               const int working_prec) {
    /* Precompute some tables of numerical values for later look-up.
     *
     * generators: List of generators
     * working_prec: Working precision
     */
    moments_t M = compute_moments(2*generators.size()+1);
    ai_t A = compute_ai(generators, M, working_prec);
    wft_t WF = compute_weightfactors(generators, A, working_prec);
    z_t Z = compute_z_sequence(A);

    return std::make_tuple(M, A, WF, Z);
}


template<int D>
nodes_t<D>
compute_nodes(const partition_t<D> P,
              const generators_t& generators,
              const int working_prec) {
    /* Compute fully symmetric quadrature nodes for given partition `P`.
     *
     * P: Partition `P`
     * generators: Table with precomputed generators
     * working_prec: Working precision
     */
    nodes_t<D> nodes;
    partition_t<D> Q;

    // Number of 0 entries in P
    int xi = nz<D>(P);
    int nsf = xi<D ? (2<<(D-xi-1)) : 1;

    permutations_t<D> permutations = Permutations<D>(P);
    for(auto it=permutations.begin(); it != permutations.end(); it++) {
        Q = *it;
        for(int v=0; v < nsf; v++) {
            node_t<D> node;
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
                node[d] = *t;
            }
            nodes.push_back(node);
        }
    }
    return nodes;
}


template<int D>
weights_t
compute_weights(const partition_t<D> P,
                const int K,
                const arb_mat_struct& weight_factors,
                const int working_prec) {
    /* Function to compute weights for partition `P`.
     *
     * P: Partition `P`
     * K: Rule order
     * weight_factors: Table with precomputed weight factors
     * working_prec: Working precision
     */
    weights_t weights;

    arb_t W, w;
    arb_init(W);
    arb_init(w);
    arb_zero(W);

    latticepoints_t<D> U = LatticePoints<D>(K-sum<D>(P));
    partition_t<D> Q;

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


template<int D>
rule_t<D>
genz_keister_construction(const int K,
                          const generators_t& generators,
                          const tables_t& tables,
                          const int working_prec) {
    /* Compute the Genz-Keister construction.
     *
     * K: Level of the quadrature rule
     * generators: Table with precomputed generators
     * weight_factors: Table with the weight factors
     * Z:
     * working_prec: Working precision
     */
    wft_t weight_factors = std::get<2>(tables);
    z_t Z = std::get<3>(tables);

    nodes_t<D> nodes;
    weights_t weights;

    // Iterate over all relevant integer partitions
    partitions_t<D> partitions = Partitions<D>(K);
    for(auto it=partitions.begin(); it != partitions.end(); it++) {
        //
        partition_t<D> P = *it;
        int s = 0;
        for(int d=0; d < D; d++) {
            s += P[d];
            s += Z[P[d]];
        }
        //
        if(s <= K) {
            // Compute nodes and weights for given partition
            nodes_t<D> p = compute_nodes<D>(P, generators, working_prec);
            weights_t w = compute_weights<D>(P, K, weight_factors, working_prec);
            for(int i=0; i < p.size(); i++) {
                nodes.push_back(p[i]);
                weights.push_back(w[0]);
            }
        }
    }

    return std::make_pair(nodes, weights);
}


bool
check_accuracy(const arb_t a,
               const long prec) {
    /* Check if a ball has a radius small enough to fit the target precision.
     *
     * a: A ball to test
     * prec: Target precision in number bits
     */
    if(mag_cmp_2exp_si(arb_radref(a), -prec) >= 0) {
        return false;
    }
    return true;
}


template<int D>
bool
check_accuracy(const nodes_t<D>& nodes,
               const weights_t& weights,
               const int target_prec) {
    /* Check if the nodes and weights are accurate enough to fit
     * the target precision.
     */
    // Check nodes
    for(auto it=nodes.begin(); it != nodes.end(); it++) {
        for(int d=0; d < D; d++) {
            if(!check_accuracy(&((*it)[d]), target_prec)) {
                return false;
            }
        }
    }
    // Check weights
    for(auto it=weights.begin(); it != weights.end(); it++) {
        if(!check_accuracy(&(*it), target_prec)) {
            return false;
        }
    }
    return true;
}


#endif
