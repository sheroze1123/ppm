#include "common.h"

#include <cstring>
#include <fftw3.h>
#include <random>

using namespace std;

const double G = 6.6748 * 10e-11;

void random_particle_initialization(const int N_p,
                                    const double pos_min, const int pos_max,
                                    const double mass_min, const int mass_max,
                                    double* restrict p_vel,
                                    double* restrict p_pos,
                                    double* restrict p_mass) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> pos_dis(pos_min, pos_max);
    uniform_real_distribution<> mass_dis(mass_min, mass_max);

    for (int i=0; i<N_p; i++) {
        p_vel[2*i]   = 0.0;
        p_vel[2*i+1] = 0.0;
        p_pos[2*i]   = pos_dis(gen);
        p_pos[2*i+1] = pos_dis(gen);
        p_mass[i]    = mass_dis(gen);
    }
}

void random_particle_initialization(const int N_p, const double L,
                                    double* restrict p_vel,
                                    double* restrict p_pos,
                                    double* restrict p_mass) {
    random_particle_initialization(N_p, 4*L/10, 6*L/10, 0.1, 3.0,
                                   p_vel, p_pos, p_mass);
}

void compute_rho(int N_p, int N, double delta_d,
                 double* restrict rho,
                 double* restrict particle_mass,
                 double* restrict particle_pos) {

    // Zero-out rho for reuse
    memset(rho, 0, N*N*sizeof(double));

    // NOTE(mjw297): this loop isn't vectorizable due to random access pattern
    int ind_x, ind_y;
    for (int i=0; i<N_p; i++) {
        ind_x =  floor(particle_pos[2*i]/delta_d);
        ind_y =  floor(particle_pos[2*i+1]/delta_d);
        rho[ind_y*N + ind_x] += particle_mass[i];
    }

    // Scale to represent density
    double scaling_factor = delta_d * delta_d;
    for (int i=0; i<N*N; i++) {
        rho[i] /= scaling_factor;
    }
}

void compute_phi_k(int N, double L, fftw_complex* rho_k) {
    double n_i, n_j, kx_i, ky_j, k_sq, scaling_factor;

    for (int j=0; j<N; j++) {

        n_j = (j >= N/2) ?  (j - N) : j;    // Negative frequency mapping in y
        ky_j = n_j * 2 * M_PI / L;          // Spatial frequency in y

        for (int i=0; i<N; i++) {
            if (i==0 && j==0) continue;

            n_i = (i >= N/2) ? (i - N) : i; // Negative frequency mapping in x
            kx_i = n_i * 2 * M_PI / L;      // Spatial frequency in x

            k_sq = kx_i * kx_i + ky_j * ky_j;
            scaling_factor = 4 * M_PI / k_sq;

            rho_k[j*N+i][0] *= (scaling_factor);
            rho_k[j*N+i][1] *= (scaling_factor);
        }
    }
}

// Given the graviational potential (phi), compute the graviational acceleration
// at each mesh point using finite difference. 
//      \vec{g} = - \triangle{\phi} 
void compute_accelerations(int N, double* restrict acc_x,
                                  double* restrict acc_y,
                                  const double* restrict phi,
                                  const double delta_d) {

    // 2 delta_d for finite difference and N*N for FFTW normalization
    double scaling_factor = 2 * delta_d * N * N;
    for (int j=0; j<N; j++) {
        for (int i=1; i<N-1; i++) {
            acc_x[j*N + i] = (-phi[j*N + i-1] +  phi[j*N + i+1]) / scaling_factor;
        }
    }

    for (int j=1; j<N-1; j++) {
        for (int i=0; i<N; i++) {
            acc_y[j*N + i] = (-phi[(j-1)*N + i] +  phi[(j+1)*N + i]) / scaling_factor;
        }
    }

    // Boundary conditions for phi are periodic
    // Left-right edges
    for (int j=0; j<N; j++) {
        acc_x[j*N] = (-phi[j*N + (N-1)] + phi[j*N + 1]) / scaling_factor;
        acc_x[j*N+(N-1)] = (-phi[j*N + N-2] + phi[j*N]) / scaling_factor;
    }

    // Top-bottom edges
    for (int i=0; i<N; i++) {
        acc_y[i] = (-phi[(N-1)*N +i] +  phi[N + i]) / scaling_factor;
        acc_y[(N-1)*N + i] = (-phi[(N-2)*N + i] + phi[i])/ scaling_factor;
    }
}

void update_particles(int N_p, int N, double delta_t, double delta_d, double L,
                      double* restrict particle_pos,
                      double* restrict particle_vel,
                      double* restrict a_x,
                      double* restrict a_y) {

    int ind_x, ind_y;
    for (int i=0; i<N_p; i++) {
        ind_x =  floor(particle_pos[2*i]/delta_d);
        ind_y =  floor(particle_pos[2*i+1]/delta_d);
        particle_vel[2*i] += a_x[ind_y*N + ind_x] * delta_t;
        particle_vel[2*i+1] += a_y[ind_y*N + ind_x] * delta_t;
        particle_pos[2*i] += particle_vel[2*i] * delta_t;
        particle_pos[2*i+1] += particle_vel[2*i+1] * delta_t;

        // Applying periodic boundary conditions
        if (particle_pos[2*i] < 0.0) particle_pos[2*i]     = fmod(particle_pos[2*i], L) + L;
        if (particle_pos[2*i] > L) particle_pos[2*i]       = fmod(particle_pos[2*i], L);
        if (particle_pos[2*i+1] < 0.0) particle_pos[2*i+1] = fmod(particle_pos[2*i+1], L) + L;
        if (particle_pos[2*i+1] > L) particle_pos[2*i+1]   = fmod(particle_pos[2*i+1],L);
    }
}
