/*  Author: R. Bourquin
    Copyright: (C) 2014 R. Bourquin
    License: GNU GPL v2 or above
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>

#include "enumerators.h"


template <int D>
void print(std::array<int, D> K) {
    std::cout << "(";
    for(int i=0; i < D-1; i++) {
        std::cout << K[i] << ", ";
    }
    std::cout << K[D-1] << ")\n";
}


template <int D>
void print(std::list<std::array<int, D> > L) {
    for(auto it = L.begin(); it != L.end(); it++) {
        auto k = *it;
        std::cout << "(";
        for(int i=0; i < D-1; i++) {
            std::cout << k[i] << ", ";
        }
        std::cout << k[D-1] << ")\n";
    }
}


int main(int argc, char* argv[]) {

    const int D = 3;

    // Test partitions
    std::cout << "Partitions\n";
    std::cout << "----------\n";

    std::list<std::array<int, D> > PL;
    PL = Partitions<D>(3);

    print<D>(PL);


    std::cout << std::endl;


    // Test lattice points
    std::cout << "Lattice Points\n";
    std::cout << "--------------\n";

    const int E = 3;
    std::list<std::array<int, E> > LPL;
    LPL = LatticePoints<E>(3);

    print<E>(LPL);


    std::cout << std::endl;


    // Test permutations
    std::cout << "Permutations\n";
    std::cout << "------------\n";

    const int F = 3;
    std::list<std::array<int, F> > L;
    std::list<std::array<int, D> > PP;
    L = Partitions<F>(6);

    for(auto Lit = L.begin(); Lit != L.end(); Lit++) {
        auto P = *Lit;

        std::cout << "[";
        for(int i=0; i < F-1; i++) {
            std::cout << P[i] << ", ";
        }
        std::cout << P[F-1] << "]\n";

        PP = Permutations<F>(P);

        print<F>(PP);
        std::cout << "Number permutations: " << PP.size() << std::endl;
        std::cout << std::endl;
    }


    return EXIT_SUCCESS;
}
