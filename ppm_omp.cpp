#include <chrono>
#include <cmath>
#include <complex.h>
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
    "  - t -- time step (0.1)\n"
    "  - s -- number of time steps (1000)\n";

void strong_scaling(int t_steps, double* particle_pos, double* particle_vel,
                    double* particle_mass, bool* particle_valid, int N_p, int N,
                    double L, double* a_x, double* a_y, double* rho, double* phi,
                    fftw_complex* rho_k, double delta_d, double delta_t) {
    FILE *fp;
    fp = fopen("strong_scaling.csv", "w+");
    double t_start, t_end, t_threadrun;
    int thread_max = 26, i;

    for (i = 1; i <= thread_max; i++) {
        omp_set_num_threads(i);
        random_particle_initialization(N_p, L,
                                       particle_vel, particle_valid, particle_pos,
                                       particle_mass);

        fftw_plan_with_nthreads(i);
        fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
        fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

        t_threadrun = 0.0;

        for (int t=1; t<t_steps; t++) {
            t_start = omp_get_wtime();

            compute_rho(N_p, N, delta_d, rho,
                        particle_mass, particle_pos, particle_valid);

            fftw_execute(rho_plan);

            compute_phi_k(N, L, rho_k);

            fftw_execute(phi_plan);

            // TODO: Scaling delta_t

            compute_accelerations(N, a_x, a_y, phi);

            update_particles(N_p, N, delta_t, delta_d, L,
                             particle_pos, particle_vel, particle_valid, a_x, a_y);


            t_end = omp_get_wtime();
            t_threadrun += (t_end - t_start);
        }
        fftw_destroy_plan(rho_plan);
        fftw_destroy_plan(phi_plan);
        fprintf(fp, "%d, %f\n", i, (t_threadrun / (double)t_steps));
    }
    fclose(fp);
}

int main(int argc, char** argv) {
    // INITIALIZATION ////////////////////////////

    int N=128, N_p = 300, t_steps=1000;
    double L = 100.0;
    double delta_t = 0.1;

    int error_code = fftw_init_threads();
    if (error_code == 0) {
        cout << "Error initializing FFTW threads. Exiting..." << endl;
        exit( EXIT_FAILURE );
    }

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
        case 's': t_steps = atoi(optarg); break;
        case 'l': L = atof(optarg);  break;
        case 't': delta_t = atof(optarg);  break;
        }
    }
    double delta_d = L/N;

    double *rho = (double*) malloc (N * N * sizeof(double));
    double *phi = (double*) malloc (N * N * sizeof(double));
    double *a_x = (double*) malloc (N * N * sizeof(double));
    double *a_y = (double*) malloc (N * N * sizeof(double));
    bool *particle_valid  = (bool*) malloc (N_p * sizeof(bool));
    double *particle_pos  = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_vel  = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_mass = (double*) malloc (N_p * sizeof(double));
    fftw_complex *rho_k   = (fftw_complex*) fftw_malloc (N * N * sizeof(fftw_complex));

    random_particle_initialization(N_p, L,
                                   particle_vel, particle_valid, particle_pos,
                                   particle_mass);

    Marshaller marshaller("particles.txt", L, N, N_p, particle_mass);
    marshaller.marshal(particle_valid, particle_pos);

    fftw_plan_with_nthreads(omp_get_max_threads());
    fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
    fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

    // TIME STEP ////////////////////////////////
    double average_time_ms = 0.0;

    for (int t=1; t<1000; t++) {
        double t_start = omp_get_wtime();

        compute_rho(N_p, N, delta_d, rho,
                    particle_mass, particle_pos, particle_valid);

        fftw_execute(rho_plan);

        compute_phi_k(N, L, rho_k);

        fftw_execute(phi_plan);

        // TODO: Scaling delta_t

        compute_accelerations(N, a_x, a_y, phi);

        update_particles(N_p, N, delta_t, delta_d, L,
                         particle_pos, particle_vel, particle_valid, a_x, a_y);

        double t_end = omp_get_wtime();
        average_time_ms += (t_end - t_start);

        marshaller.marshal(particle_valid, particle_pos);
    }

    cout << "Average time per step in milliseconds: " << average_time_ms/t_steps << endl;

    // strong_scaling(t_steps, particle_pos, particle_vel, particle_mass, particle_valid, N_p, N,
            // L, a_x, a_y, rho, phi, rho_k, delta_d, delta_t);

    fftw_destroy_plan(rho_plan);
    fftw_destroy_plan(phi_plan);
    free(particle_pos);
    free(particle_vel);
    free(particle_mass);
    free(rho);
    free(a_x);
    free(a_y);
    free(particle_valid);
    fftw_free(rho_k);
    fftw_cleanup_threads();
}
