#ifndef COMMON_H
#define COMMON_H

#include <fftw3.h>

extern const double G; //TODO: Scale everything by G

// Performs the following initialization.
//   1. `p_vel` is set to all 0
//   2. `p_pos` is set uniformally at random between `pos_min` and `pos_max`
//   3. `p_mass` is set uniformally at random between `mass_min` and `mass_max`
void random_particle_initialization(const int N_p,
                                    const double pos_min, const int pos_max,
                                    const double mass_min, const int mass_max,
                                    double* restrict p_vel,
                                    double* restrict p_pos,
                                    double* restrict p_mass);

// An overload of the above function that chooses a sensible value for
// `pos_min`, `pos_max`, `mass_min`, and `mass_max`.
void random_particle_initialization(const int N_p, const double L,
                                    double* restrict p_vel,
                                    double* restrict p_pos,
                                    double* restrict p_mass);

// Compute the mass density at grid points (rho) using Nearest-Grid Point method
void compute_rho(int N_p, int N, double delta_d,
                 double* restrict rho,
                 double* restrict particle_mass,
                 double* restrict particle_pos);

// Compute the frequency domain representation of gravitational potential
// using the frequency domain representation of mass density.
void compute_phi_k(int N, double L, fftw_complex* rho_k);

// Given the graviational potential (phi), compute the graviational acceleration
// at each mesh point using finite difference.
void compute_accelerations(int N, double* restrict acc_x,
                                  double* restrict acc_y,
                                  const double* restrict phi,
                                  const double delta_d);

// Perform a time step update of particle positions and velocities.
void update_particles(int N_p, int N, double delta_t, double delta_d, double L,
                      double* restrict particle_pos,
                      double* restrict particle_vel,
                      double* restrict a_x,
                      double* restrict a_y);

#endif // COMMON_H
