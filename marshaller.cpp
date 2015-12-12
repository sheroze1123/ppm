#include "marshaller.h"

#include <cstdio>
#include <fstream>
#include <iostream>

// This is the marshall file format:
//
//   n                      // the number of times marshal was called
//   L                      // the side length of the simulation
//   N                      // the number of grid points
//   Np                     // the number of particles
//   mass_1                 // the particle masses
//   ...
//   mass_Np
//   x_1 y_1 vx_1 vy_1      // the particle positions and velocities
//   ...
//   x_Np y_Np vx_Np vy_Np
//   x_1 y_1 vx_1 vy_1
//   ...
//   x_Np y_Np vx_Np vy_Np

Marshaller::Marshaller(const std::string& filename, const double L, const int N,
                       const int N_p, const double *mass)
    : f_(filename), N_p_(N_p), n_(0)
{
    // The first thing in the marshal file is the number of times marshal is
    // called. Here, we preserve a bit of space at the beginning of the file by
    // writing 1000000000. When the marshaller is destructed, this is rewritten
    // with the actuall value of n.
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
    char n[11] = {0};
    sprintf(n, "%010d", n_);
    f_.seekp(0);
    f_ << n << std::endl;
}

void Marshaller::marshal(const bool* valid, const double *positions) {
    // TODO: fix
    for (int i = 0; i < N_p_; ++i) {
        f_ << positions[2*i    ] << " "
           << positions[2*i + 1];
        f_.flush();
    }
    ++n_;
}
