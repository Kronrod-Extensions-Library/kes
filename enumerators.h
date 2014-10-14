/*  Author: R. Bourquin
 *  Copyright: (C) 2014 R. Bourquin
 *  License: GNU GPL v2 or above
 *
 *  A library of helper functions to search for
 *  Kronrod extensions of Gauss quadrature rules.
 */

#ifndef __HH__enumerators
#define __HH__enumerators

#include <array>
#include <list>


template<int D> using partition_t = std::array<int, D>;
template<int D> using partitions_t = std::list<partition_t<D>>;
template<int D> using point_t = std::array<int, D>;
template<int D> using latticepoints_t = std::list<point_t<D>>;
template<int D> using permutation_t = std::array<int, D>;
template<int D> using permutations_t = std::list<permutation_t<D>>;


template<int D>
int sum(const point_t<D> Z) {
    /* Sum of the components of the point Z.
     */
    int s = 0;
    for(int i=0; i<D; i++) {
        s += Z[i];
    }
    return s;
}


template<int D>
int nz(const point_t<D> Z) {
    /* Number of zero components of the point Z.
     */
    int s = 0;
    for(int i=0; i<D; i++) {
        if(Z[i] == 0) {
            s++;
        }
    }
    return s;
}


template<int D>
int nnz(const point_t<D> Z) {
    /* Number of non-zero components of the point Z.
     */
    return D - nz<D>(Z);
}


template<int D>
partitions_t<D>
Partitions(const int K) {
    /*
     * Enumerate integer partitions in anti-lexocographic
     * order for integers up to some limit K. All partitions
     * have exactly D parts, some may be zero.
     *
     * K: Level
     */
    partitions_t<D> partitions;
    partition_t<D> P;

    P.fill(0);
    partitions.push_back(P);

    while(sum<D>(P) <= K) {
        int p0 = P[0];
        bool broke = false;
        for(int i=1; i < D; i++) {
            p0 += P[i];
            if(P[0] <= P[i] + 1) {
                P[i] = 0;
            } else {
                P[0] = p0 - i * (P[i] + 1);
                for(int j=1; j <= i; j++) {
                    P[j] = P[i] + 1;
                }
                partitions.push_back(P);
                broke = true;
                break;
            }
        }
        if(!broke) {
            P[0] = p0 + 1;
            if(sum<D>(P) <= K) {
                partitions.push_back(P);
            }
        }
    }

    return partitions;
}


template<int D>
latticepoints_t<D>
LatticePoints(const int N) {
    /* This method enumerates all lattice points of a lattice
     * :math:`\Lambda \subset \mathbb{N}^D` in :math:`D` dimensions
     * having fixed :math:`l_1` norm :math:`N`.
     *
     * N: The :math:`l_1` norm of the lattice points.
     */
    latticepoints_t<D> L;
    point_t<D> k;

    for(int n=0; n <= N; n++) {
        k.fill(0);
        k[0] = n;
        L.push_back(k);
        int c = 1;
        while(k[D-1] < n) {
            if(c == D) {
                for(int i = c-1; i > 0; i--) {
                    c = i;
                    if(k[i-1] != 0) {
                        break;
                    }
                }
            }
            k[c-1]--;
            c++;
            k[c-1] = n;
            for(int t=0; t < c-1; t++) {
                k[c-1] -= k[t];
            }
            if(c < D) {
                for(int t=c; t < D; t++) {
                    k[t] = 0;
                }
            }
            L.push_back(k);
        }
    }

    return L;
}


template<int D>
permutations_t<D>
Permutations(const permutation_t<D> permutation) {
    /* Enumerate all permutations in anti-lexicographical
     * order follwing the given permutation `P`.
     *
     * P: A permutation
     */
    permutations_t<D> permutations;
    permutation_t<D> P = permutation;

    permutations.push_back(P);

    bool broke = true;
    while(broke) {
        broke = false;
        for(int i=1; i < D; i++) {
            int pi = P[i];
            if(P[i-1] > pi) {
                int I = i;
                if(i > 1) {
                    int J = I;
                    int u = I / 2;
                    for(int j=0; j < u; j++) {
                        int pj = P[j];
                        if(pj <= pi) {
                            I--;
                        }
                        P[j] = P[i-j-1];
                        P[i-j-1] = pj;
                        if(P[j] > pi) {
                            J = j + 1;
                        }
                    }
                    if(P[I-1] <= pi) {
                        I = J;
                    }
                }
                P[i] = P[I-1];
                P[I-1] = pi;
                permutations.push_back(P);
                broke = true;
                break;
            }
        }
    }

    return permutations;
}


#endif
