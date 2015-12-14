#include "marshaller.h"

#include <cstdio>
#include <fstream>
#include <iostream>

// This is the marshall file format:
//
//   n         // the number of time steps
//   L         // specifies an L by L simulation space
//   N         // specifies an N by N grid
//   Np        // the number of particles
//   mass_1    // the particle masses
//   ...
//   mass_Np
//   x_1 y_1   // the particle positions, or "-1 -1" if the
//   x_2 y_2   // particle is out of bounds
//   ...
//   x_Np y_Np
//   x_1 y_1
//   ...
//   x_Np y_Np

Marshaller::Marshaller(const std::string& filename, const double L, const int N,
                       const int N_p, const double *mass)
    : f_(filename), N_p_(N_p), n_(0)
{
    // The first thing in the marshal file is the number of time steps. Here,
    // we preserve a bit of space at the beginning of the file by writing
    // 1000000000. When the marshaller is destructed, this is rewritten with
    // the actuall value of n_.
    f_ << "1000000000" << std::endl;
    f_ << L << std::endl;
    f_ << N << std::endl;
    f_ << N_p << std::endl;
    for (int i = 0; i < N_p_; ++i) {
        f_ << mass[i] << std::endl;
    }
    f_.flush();
}

Marshaller::~Marshaller() {
    // Here, we overwrite the beginning of the file with the actual number of
    // time steps. 11 is the length of the null terminated string 1000000000.
    // The number is written with leading padding zeros to make it fit exactly.
    char n[11] = {0};
    sprintf(n, "%010d", n_);
    f_.seekp(0);
    f_ << n << std::endl;
}

void Marshaller::marshal(const double *positions) {
    for (int i = 0; i < N_p_; ++i) {
        f_ << positions[2*i] << " " << positions[2*i + 1] << std::endl;
        f_.flush();
    }
    ++n_;
}
