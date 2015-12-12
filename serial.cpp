#include <getopt.h>
#include <cmath>
#include <complex.h>
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "marshaller.h"

using namespace std;

const double G = 6.6748 * 10e-11; //TODO: Scale everything by G

const char* usage =
    "serial -- Serial N-body simulation using a particle-mesh method\n"
    "Flags:\n"
    "  - n -- number of grid points of the mesh\n"
    "  - p -- number of particles\n"
    "  - l -- side length of simulation\n"
    "  - t -- time step (0.1)\n";

int main(int argc, char** argv) {
    // INITIALIZATION ////////////////////////////

    int N=128, N_p = 300;
    double L = 100.0;
    double delta_t = 0.1;

    // Option processing
    extern char* optarg;
    const char* optstring = "hp:n:l:t:";
    int c;
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch (c) {
        case 'h':
            fprintf(stderr, "%s", usage);
            return -1;
        case 'n': N = atoi(optarg); break;
        case 'p': N_p = atoi(optarg); break;
        case 'l': L = atof(optarg);  break;
        case 't': delta_t = atof(optarg);  break;
        }
    }
    double delta_d = L/N;

    double *rho = (double*) malloc (N * N * sizeof(double));
    double *phi = (double*) malloc (N * N * sizeof(double));
    double *a_x = (double*) malloc (N * N * sizeof(double));
    double *a_y = (double*) malloc (N * N * sizeof(double));
    fftw_complex *rho_k = (fftw_complex*) fftw_malloc (N * N * sizeof(fftw_complex));
    double *particle_pos = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_vel = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_mass = (double*) malloc (N_p * sizeof(double));

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> pos_dis(40.0, 60.0);
    uniform_real_distribution<> mass_dis(0.1, 3.0);
    for (int i=0; i<N_p; i++) {
        particle_pos[2*i]   = pos_dis(gen);
        particle_pos[2*i+1] = pos_dis(gen);
        particle_vel[2*i]   = 0.0;
        particle_vel[2*i+1] = 0.0;
        particle_mass[i]    = mass_dis(gen);
    }

    // Marshaller marshaller("particles.txt", L, N, N_p, particle_mass);
    // marshaller.marshal(NULL [>TODO<], particle_pos);

    fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_ESTIMATE);
    fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_ESTIMATE);

    // TIME STEP ////////////////////////////////
    double time = 0.0;
    for (int t=1; t<1000; t++) {

        time = t*delta_t;

        for (int i=0; i<N*N; i++) {
            rho[i] = 0.0;
        }
        int ind_x, ind_y;
        for (int i=0; i<N_p; i++) {
            ind_x =  floor(particle_pos[2*i]/delta_d);
            ind_y =  floor(particle_pos[2*i+1]/delta_d);
            rho[ind_y*N + ind_x] += particle_mass[i];
        }

        for (int i=0; i<N*N; i++) {
            rho[i] /= (delta_d * delta_d);
        }

        fftw_execute(rho_plan);

        double n_i, n_j, kx_i, ky_j, k_sq, mult, y_mult, a;
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
        fftw_execute(phi_plan);
        // double max = 0, E;
        // for (int i=0; i<N_p; i++) {
            // ind_x =  floor(particle_pos[2*i]/delta_d);
            // ind_y =  floor(particle_pos[2*i+1]/delta_d);
            // E = E_x[ind_y*N + ind_x];
            // if (E > max) {
                // max = E;
            // }
            // E = E_y[ind_y*N + ind_x];
            // if (E > max) {
                // max = E;
            // }
        // }

        // delta_t = sqrt(L/max);

        // TODO: Edge cases
        for (int j=1; j<N-1; j++) {
            for (int i=1; i<N-1; i++) {
                a_x[j*N + i] = (phi[j*N + i-1] -  phi[j*N + i+1])/(2*N*N);
                a_y[j*N + i] = (phi[(j-1)*N + i] -  phi[(j+1)*N + i])/(2*N*N);
            }
        }

        for (int i=0; i<N_p; i++) {
            ind_x =  floor(particle_pos[2*i]/delta_d);
            ind_y =  floor(particle_pos[2*i+1]/delta_d);
            particle_vel[2*i] -= a_x[ind_y*N + ind_x] * delta_t;
            particle_vel[2*i+1] -= a_y[ind_y*N + ind_x] * delta_t;
            particle_pos[2*i] += particle_vel[2*i] * delta_t;
            particle_pos[2*i+1] += particle_vel[2*i+1] * delta_t;

            // TODO: Particles out of bounds handling
            if (particle_pos[2*i] < 0.0) particle_pos[2*i] += L;
            if (particle_pos[2*i] > L) particle_pos[2*i] -= L;
            if (particle_pos[2*i+1] < 0.0) particle_pos[2*i+1] += L;
            if (particle_pos[2*i+1] > L) particle_pos[2*i+1] -= L;
        }
        // marshaller.marshal(NULL /* TODO */, particle_pos);
    }


    fftw_destroy_plan(rho_plan);
    fftw_destroy_plan(phi_plan);
    free(particle_pos);
    free(particle_vel);
    free(particle_mass);
    free(rho);
    free(a_x);
    free(a_y);
    fftw_free(rho_k);
}
