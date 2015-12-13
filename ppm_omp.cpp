#include <omp.h>
#include <getopt.h>
#include <cmath>
#include <complex.h>
#include <fftw3.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "marshaller.h"

using namespace std;

const double G = 6.6748 * 10e-11; //TODO: Scale everything by G

// Initialized the given particles to with random positions and masses but with zero initial velocity.
void random_particle_initialization(int N_p, double* p_pos, double* p_vel, double* p_mass, bool* p_valid) {

    random_device rd;
    mt19937 gen(rd());
    double pos_min = 30.0, pos_max = 70.0, mass_min = 0.1, mass_max = 3.0;
    uniform_real_distribution<> pos_dis(pos_min, pos_max);
    uniform_real_distribution<> mass_dis(mass_min, mass_max);

    for (int i=0; i<N_p; i++) {
        p_valid[i] = true;
        p_pos[2*i]   = pos_dis(gen);
        p_pos[2*i+1] = pos_dis(gen);
        p_vel[2*i]   = 0.0;
        p_vel[2*i+1] = 0.0;
        p_mass[i]    = mass_dis(gen);
    }
}

// Given the graviational potential (phi), compute the graviational acceleration
// at each mesh point using finite difference.
void compute_accelerations(int N, double* acc_x, double* acc_y, double* phi) {
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

// Perform a time step update of particle positions and velocities.
void update_particles(int N_p, int N, double* particle_pos, double* particle_vel, bool* particle_valid,
        double* a_x, double* a_y, double delta_t, double delta_d, double L) {

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

// Compute the mass density at grid points (rho) using Nearest-Grid Point method
void compute_rho(int N_p, int N, double *rho, double* particle_mass, double* particle_pos,
        bool* particle_valid, double delta_d) {

    for (int i=0; i<N*N; i++) {
        rho[i] = 0.0;
    }
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

// Compute the frequency domain representation of gravitational potential
// using the frequency domain representation of mass density.
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

const char* usage =
    "serial -- Serial N-body simulation using a particle-mesh method\n"
    "Flags:\n"
    "  - n -- number of grid points of the mesh\n"
    "  - p -- number of particles\n"
    "  - l -- side length of simulation\n"
    "  - t -- time step (0.1)\n"
    "  - s -- number of time steps (1000)\n";

void strong_scaling(int t_steps, double* particle_pos, double* particle_vel, double* particle_mass,
        bool* particle_valid, int N_p, int N, double L, double* a_x, double* a_y, double* rho, double* phi,
        fftw_complex* rho_k, double delta_d, double delta_t) {

    FILE *fp;
    fp = fopen("strong_scaling.csv", "w+");
    double t_start, t_end, t_threadrun;
    int thread_max = 26, i;

    for (i = 1; i <= thread_max; i++) {
        omp_set_num_threads(i);
        random_particle_initialization(N_p, particle_pos, particle_vel, particle_mass, particle_valid);

        fftw_plan_with_nthreads(i);
        fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
        fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

        t_threadrun = 0.0;

        for (int t=1; t<t_steps; t++) {
            t_start = omp_get_wtime();

            compute_rho(N_p, N, rho, particle_mass, particle_pos, particle_valid, delta_d);

            fftw_execute(rho_plan);

            compute_phi_k(N, L, rho_k);

            fftw_execute(phi_plan);

            // TODO: Scaling delta_t

            compute_accelerations(N, a_x, a_y, phi);

            update_particles(N_p, N, particle_pos, particle_vel, particle_valid, a_x, a_y, delta_t, delta_d, L);


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

    random_particle_initialization(N_p, particle_pos, particle_vel, particle_mass, particle_valid);
    Marshaller marshaller("particles.txt", L, N, N_p, particle_mass);
    marshaller.marshal(particle_valid, particle_pos);

    fftw_plan_with_nthreads(omp_get_max_threads());
    fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
    fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

    // TIME STEP ////////////////////////////////
    double average_time_ms = 0.0;

    for (int t=1; t<1000; t++) {
        double t_start = omp_get_wtime();

        compute_rho(N_p, N, rho, particle_mass, particle_pos, particle_valid, delta_d);

        fftw_execute(rho_plan);

        compute_phi_k(N, L, rho_k);

        fftw_execute(phi_plan);

        // TODO: Scaling delta_t

        compute_accelerations(N, a_x, a_y, phi);

        update_particles(N_p, N, particle_pos, particle_vel, particle_valid, a_x, a_y, delta_t, delta_d, L);

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
