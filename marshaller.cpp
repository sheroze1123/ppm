#include "marshaller.h"

#include <fstream>
#include <iostream>

Marshaller::Marshaller(const std::string& filename, const double L, const int N,
                       const int N_p, const double *mass)
    : f_(filename), N_p_(N_p)
{
    f_ << L << std::endl;
    f_ << N << std::endl;
    f_ << N_p << std::endl;
    for (int i = 0; i < N_p_; ++i) {
        f_ << mass[i] << std::endl;
    }
    f_.flush();
}

void Marshaller::marshal(const double *position, const double *velocity) {
    for (int i = 0; i < N_p_; ++i) {
        f_ << position[2*i    ] << " "
           << position[2*i + 1] << " "
           << velocity[2*i    ] << " "
           << velocity[2*i + 1] << std::endl;
        f_.flush();
    }
}
