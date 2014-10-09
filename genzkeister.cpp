/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <array>
#include <list>
#include <iostream>

//#include "libkes2.h"

template <int D>
int sum(std::array<int, D> Z) {
    int s = 0;
    for(int i=0; i<D; i++) {
	s += Z[i];
    }
    return s;
}


template <int D>
std::list<std::array<int, D> >
Partitions(int K) {
    /*
     * Enumerate integer partitions in anti-lexocographic
     * order for integers up to some limit K. All partitions
     * have exactly D parts, some may be zero.
     *
     * :param D: Dimension
     * :param K: Level
     */
    std::list<std::array<int, D> > partitions;
    std::array<int, D> P;

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


template <int D>
std::list<std::array<int, D> >
LatticePoints(int N) {
    /* This method enumerates all lattice points of a lattice
     * :math:`\Lambda \subset \mathbb{N}^D` in :math:`D` dimensions
     * having fixed :math:`l_1` norm :math:`N`.
     *
     * :param D: The dimension :math:`D` of the lattice.
     * :param N: The :math:`l_1` norm of the lattice points.
     */
    std::list<std::array<int, D> > L;
    std::array<int, D> k;

    for(int n=0; n <= N; n++) {
	k.fill(0);
	k[0] = n;
	L.push_back(k);

	int c = 1;
	while(k[D-1] < n) {
	    if(c == D) {
		for(int i = c-1; i >= 1; i--) {
		    c = i;
		    if(k[i-1] != 0) {
			break;
		    }
		}
	    }
	    k[c-1] -= 1;
	    c += 1;

	    int xi = 0;
	    for(int t=0; t < c-1; t++) {
		xi += k[t];
	    }
	    k[c-1] = n - xi;

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



template <int D>
std::list<std::array<int, D> >
Permutations(std::array<int, D> P) {
    /* Enumerate all permutations in anti-lexicographical
     * order follwing the given permutation `P`.
     *
     * :param P: A permutation
     */
    std::list<std::array<int, D> > permutations;

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
		    for(int j=0; j < I/2; j++) {
			int pj = P[j];
			if(pj <= pi) {
			    I = I - 1;
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




int main(int argc, char* argv[]) {
    /* int i, j; */
    /* int deg; */
    /* //fmpq_poly_t Pn; */
    /* char *strf; */
    /* //  acb_ptr nodes; */
    /* //acb_ptr weights; */
    /* int working_prec; */
    /* int target_prec; */
    /* int nrprintdigits; */
    /* int loglevel; */

    if(argc <= 1) {
        printf("Compute Genz-Keister quadrature rule\n");
        printf("Syntax: gk [-dc D] [-dp D] [-l L] n\n");
        printf("Options:\n");
        printf("        -dc  Compute nodes and weights up to this number of decimal digits\n");
        printf("        -dp  Print this number of decimal digits\n");
        printf("        -l   Set the log level\n");
        return EXIT_FAILURE;
    }


    // Test lattice point code
    const int D = 4;
    std::list<std::array<int, D> > LL;
    LL = Partitions<D>(5);

    for(auto Xit = LL.begin(); Xit != LL.end(); Xit++) {
	auto X = *Xit;


	std::cout << "[";
	for(int i=0; i < D-1; i++) {
	    std::cout << X[i] << ", ";
	}
	std::cout << X[D-1] << "]\n";

	std::list<std::array<int, D> > PP;
	PP = Permutations<D>(X);

	for (auto it = PP.begin(); it != PP.end(); it++) {
	    auto k = *it;
	    std::cout << "(";
	    for(int i=0; i < D-1; i++) {
		std::cout << k[i] << ", ";
	    }
	    std::cout << k[D-1] << ")\n";
	}
	std::cout << "\n";

    }



    return EXIT_SUCCESS;
}
