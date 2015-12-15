#include <chrono>
#include <cmath>
#include <complex.h>
#include <cstring>
#include <fftw3.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#include "common.h"
#include "marshaller.h"

using namespace std;

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
    int T = 1000;

    // Option processing
    extern char* optarg;
    const char* optstring = "hp:n:l:t:s:";
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
        case 's': T = atof(optarg);  break;
        }
    }
    double delta_d = L/N;

    double *rho = (double*) _mm_malloc (N * N * sizeof(double), 64);
    double *phi = (double*) _mm_malloc (N * N * sizeof(double), 64);
    double *a_x = (double*) _mm_malloc (N * N * sizeof(double), 64);
    double *a_y = (double*) _mm_malloc (N * N * sizeof(double), 64);
    double *particle_pos  = (double*) _mm_malloc (N_p * 2 * sizeof(double), 64);
    double *particle_vel  = (double*) _mm_malloc (N_p * 2 * sizeof(double), 64);
    double *particle_mass = (double*) _mm_malloc (N_p * sizeof(double), 64);
    fftw_complex *rho_k   = (fftw_complex*) fftw_malloc (N * N * sizeof(fftw_complex));

    random_particle_initialization(N_p, L,
                                   particle_vel, particle_pos,
                                   particle_mass);

#ifdef MARSHAL
    Marshaller marshaller("particles.txt", L, N, N_p, particle_mass);
    marshaller.marshal(particle_pos);
#endif

    fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
    fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

    // TIME STEP ////////////////////////////////
    double t_start = omp_get_wtime();

    for (int t=1; t<T; t++) {

        compute_rho(N_p, N, delta_d, rho,
                    particle_mass, particle_pos);

        fftw_execute(rho_plan);

        compute_phi_k(N, L, rho_k);

        fftw_execute(phi_plan);

        // TODO: Scaling delta_t

        compute_accelerations(N, a_x, a_y, phi, delta_d);

        update_particles(N_p, N, delta_t, delta_d, L,
                         particle_pos, particle_vel, a_x, a_y);

#ifdef MARSHAL
        marshaller.marshal(particle_pos);
#endif
    }

    double t_end = omp_get_wtime();
    double total_time_ms = t_end - t_start;
    double average_time_ms = total_time_ms / T;

    cout << "serial "
         << L << " "
         << N << " "
         << N_p << " "
         << delta_t << " "
         << T << " "
         << 1 << " "
         << total_time_ms << " "
         << average_time_ms << " "
         << endl;

    fftw_destroy_plan(rho_plan);
    fftw_destroy_plan(phi_plan);
    _mm_free(particle_pos);
    _mm_free(particle_vel);
    _mm_free(particle_mass);
    _mm_free(rho);
    _mm_free(a_x);
    _mm_free(a_y);
    fftw_free(rho_k);
}
