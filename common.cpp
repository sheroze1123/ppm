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
                                    bool*   restrict p_valid,
                                    double* restrict p_pos,
                                    double* restrict p_mass) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> pos_dis(pos_min, pos_max);
    uniform_real_distribution<> mass_dis(mass_min, mass_max);

    for (int i=0; i<N_p; i++) {
        p_vel[2*i]   = 0.0;
        p_vel[2*i+1] = 0.0;
        p_valid[i]   = true;
        p_pos[2*i]   = pos_dis(gen);
        p_pos[2*i+1] = pos_dis(gen);
        p_mass[i]    = mass_dis(gen);
    }
}

void random_particle_initialization(const int N_p, const double L,
                                    double* restrict p_vel,
                                    bool*   restrict p_valid,
                                    double* restrict p_pos,
                                    double* restrict p_mass) {
    random_particle_initialization(N_p, 4*L/10, 6*L/10, 0.1, 3.0,
                                   p_vel, p_valid, p_pos, p_mass);
}

void compute_rho(int N_p, int N, double delta_d,
                 double* restrict rho,
                 double* restrict particle_mass,
                 double* restrict particle_pos,
                 bool* restrict particle_valid) {
    memset(rho, 0, N*N*sizeof(double));

    // NOTE(mjw297): this loop isn't vectorizable due to random access pattern
    int ind_x, ind_y;
    for (int i=0; i<N_p; i++) {
        if (particle_valid[i]) {
            ind_x =  floor(particle_pos[2*i]/delta_d);
            ind_y =  floor(particle_pos[2*i+1]/delta_d);
            rho[ind_y*N + ind_x] += particle_mass[i];
        }
    }

    for (int i=0; i<N*N; i++) {
        rho[i] /= (delta_d * delta_d);
    }
}

void compute_phi_k(int N, double L, fftw_complex* rho_k) {
    double n_i, n_j, kx_i, ky_j, k_sq, mult;
    for (int j=0; j<N; j++) {

        if (j >= N/2)  n_j = (j - N); // Negative frequency mapping in y
        else n_j = j;

        ky_j = n_j * 2 * M_PI / L; // Spatial frequency in y

        for (int i=0; i<N; i++) {
            if (i==0 && j==0) continue;

            if (i >= N/2)  n_i = (i - N); // Negative frequency mapping in x
            else n_i = i;

            kx_i = n_i * 2 * M_PI / L; // Spatial frequency in x
            k_sq = kx_i * kx_i + ky_j * ky_j;
            mult = 4 * M_PI / k_sq;

            rho_k[j*N+i][0] *= (mult);
            rho_k[j*N+i][1] *= (mult);
        }
    }
}

void compute_accelerations(int N, double* restrict acc_x,
                                  double* restrict acc_y,
                                  const double* restrict phi) {
    for (int j=1; j<N-1; j++) {
        for (int i=1; i<N-1; i++) {
            acc_x[j*N + i] = (phi[j*N + i-1] -  phi[j*N + i+1])/(2*N*N);
            acc_y[j*N + i] = (phi[(j-1)*N + i] -  phi[(j+1)*N + i])/(2*N*N);
        }
    }

    // Left-right edges
    for (int j=0; j<N; j++) {
        acc_x[j*N] = -phi[j*N + 1]/(2*N*N);
        acc_x[j*N+(N-1)] = phi[j*N + N-2]/(2*N*N);
    }

    // Top-bottom edges
    for (int i=0; i<N; i++) {
        acc_x[i] = -phi[N + i]/(2*N*N);
        acc_x[(N-1)*N + i] = phi[(N-2)*N + i]/(2*N*N);
    }
}

void update_particles(int N_p, int N, double delta_t, double delta_d, double L,
                      double* restrict particle_pos,
                      double* restrict particle_vel,
                      bool*  restrict particle_valid,
                      double* restrict a_x,
                      double* restrict a_y) {
    int ind_x, ind_y;
    for (int i=0; i<N_p; i++) {
        if (particle_valid[i]) {
            ind_x =  floor(particle_pos[2*i]/delta_d);
            ind_y =  floor(particle_pos[2*i+1]/delta_d);
            particle_vel[2*i] -= a_x[ind_y*N + ind_x] * delta_t;
            particle_vel[2*i+1] -= a_y[ind_y*N + ind_x] * delta_t;
            particle_pos[2*i] += particle_vel[2*i] * delta_t;
            particle_pos[2*i+1] += particle_vel[2*i+1] * delta_t;

            if ((particle_pos[2*i] < 0.0) || (particle_pos[2*i] > L)
                || (particle_pos[2*i+1] < 0.0) ||  (particle_pos[2*i+1] > L)) {

                particle_valid[i] = false;
            }
        }
    }
}
