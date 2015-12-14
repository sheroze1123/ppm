#include <getopt.h>
#include <cmath>
#include <complex.h>
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include <string>
#include <vector>
#include <cstring>

#include "marshaller.h"

using namespace std;

const double G_scaled = 39.5; // Time in years, distances in AU, mass in solar masses

void debug_initialization(int N_p, double* p_pos, double* p_vel, double* p_mass, double L) {
    p_pos[0] = L/2;
    p_pos[1] = L/2;
    p_vel[0] = 0.0;
    p_vel[1] = 0.0;
    p_mass[0] = 1000.0;

    p_pos[2] = 2*L/3;
    p_pos[3] = L/2;
    p_vel[2] = 0.0;
    p_vel[3] = sqrt(G_scaled * 1000 * (6 / L));
    p_mass[1] = 0.001;
}

// Distance from Sun to Pluto at its most distant, 1 AU
// Radius of Earth 6380km -> 4e-5 AU
// Distance from Earth to Moon -> 1e-3 AU
// Ideal L -> 12.0 AU (Covers up to Saturn)
// Ideal N -> 1000 due to computational cost
// Ideal delta_d -> 12/1000AU -> 1.8 million km
void solar_sys_initialization(int N_p, double* p_pos, double* p_vel, double* p_mass, double L) {

    // Random number generator
    random_device rd;
    mt19937 gen(rd());

    // Mass distribution is a Poisson distribution with a mean of 1/10 of Earth
    // double mass_mean = 3.0e-7;
    double mass_mean = 1.0;
    poisson_distribution<> mass_dis(mass_mean);

    // Mean velocity of Earth 6.324 AU/year
    double vel_mean = 10.0;
    double vel_sigma = 2;
    poisson_distribution<> vel_dis(vel_mean);

    // Initial position distribution is a Gaussian centered at (L/2, L/2)
    double r_mean = L/2;
    double r_sigma = 2; //TODO: Figure out best distribution
    double x, y, r;
    normal_distribution<> gauss(r_mean, r_sigma);

    for (int i=0; i<N_p; i++) {

        p_mass[i] = mass_dis(gen);

        x = -1.0;
        y = -1.0;
        while (x < 0.0 || x > L) {
            x = gauss(gen);
        }
        while (y < 0.0 || y > L) {
            y = gauss(gen);
        }
        p_pos[2*i]   = x;
        p_pos[2*i+1] = y;
        r = sqrt(x*x + y*y);

        // s -> y value centered about L/2
        // t -> x value centered about L/2
        double s = y - L/2;
        double t = x - L/2;
        double denom = sqrt(s*s + t*t);

        // double v = vel_dis(gen);
        double v = sqrt(G_scaled * 1000 / r); 

        p_vel[2*i]   = -v * s / denom;
        p_vel[2*i+1] =  v * t / denom;
    }

    p_mass[0] = 1000.0;
    p_pos[0] = L/2;
    p_pos[1] = L/2;
    p_vel[0] = 0;
    p_vel[1] = 0;
}

// Initialized the given particles to with random positions and masses but with zero initial velocity.
void random_particle_initialization(int N_p, double* p_pos, double* p_vel, double* p_mass, double L) {

    random_device rd;
    mt19937 gen(rd());
    double pos_min = L * 0.4, pos_max = L * 0.6, mass_min = 0.1, mass_max = 3.0;
    uniform_real_distribution<> pos_dis(pos_min, pos_max);
    uniform_real_distribution<> mass_dis(mass_min, mass_max);

    for (int i=0; i<N_p; i++) {
        p_pos[2*i]   = pos_dis(gen);
        p_pos[2*i+1] = pos_dis(gen);
        p_vel[2*i]   = 0.0;
        p_vel[2*i+1] = 0.0;
        p_mass[i]    = mass_dis(gen);
    }
}

// Given the graviational potential (phi), compute the graviational acceleration
// at each mesh point using finite difference. 
//      \vec{g} = - \triangle{\phi} 
void compute_accelerations(int N, double* acc_x, double* acc_y, double* phi, double delta_d) {

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

// Perform a time step update of particle positions and velocities.
void update_particles(int N_p, int N, double* particle_pos, double* particle_vel, 
        double* a_x, double* a_y, double delta_t, double delta_d, double L) {

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

// Compute the mass density at grid points (rho) using Nearest-Grid Point method
void compute_rho(int N_p, int N, double *rho, double* particle_mass, double* particle_pos, double delta_d) {

    // Zero-out rho for reuse
    memset(rho, 0, N*N*sizeof(double));

    // Assign cell centers mass using NGP scheme
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

// Compute the frequency domain representation of gravitational potential
// using the frequency domain representation of mass density.
void compute_phi_k(int N, double L, fftw_complex* rho_k, double G) {

    double n_i, n_j, kx_i, ky_j, k_sq, scaling_factor;

    for (int j=0; j<N; j++) {

        n_j = (j >= N/2) ?  (j - N) : j;    // Negative frequency mapping in y
        ky_j = n_j * 2 * M_PI / L;          // Spatial frequency in y

        for (int i=0; i<N; i++) {
            if (i==0 && j==0) continue;

            n_i = (i >= N/2) ? (i - N) : i; // Negative frequency mapping in x
            kx_i = n_i * 2 * M_PI / L;      // Spatial frequency in x

            k_sq = kx_i * kx_i + ky_j * ky_j;
            scaling_factor = 4 * M_PI * G / k_sq;

            rho_k[j*N+i][0] *= (scaling_factor);
            rho_k[j*N+i][1] *= (scaling_factor);
        }
    }
}

const char* usage =
    "serial -- Serial N-body simulation using a particle-mesh method\n"
    "Flags:\n"
    "  - n -- number of grid points of the mesh (1024)\n"
    "  - p -- number of particles (100)\n"
    "  - l -- side length of simulation (60.0AU)\n"
    "  - t -- time step (0.0001)\n"
    "  - d -- debug mode (random particles initialization)\n"
    "  - s -- total time steps (100)\n";


int main(int argc, char** argv) {

    // INITIALIZATION ////////////////////////////

    bool debug = false;
    int N=1024, N_p = 3000, T=400;
    double L = 60.0;
    double delta_t = 0.0001;
    double G = G_scaled;

    // Option processing
    extern char* optarg;
    const char* optstring = "hp:n:l:t:d:s:";
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
        case 'd': debug = atoi(optarg); break;
        case 's': T = atoi(optarg); break;
        }
    }

    if(debug) {
        // G=1.0;
        // N=128;
        // N_p = 100;
        // L = 100.0;
        // delta_t = 0.1;

        L = 60.0;
        delta_t = 0.0005;
        G = G_scaled;
        N=1024;
        N_p = 2; 
        T=300;
    }

    double delta_d = L/N;

    double *rho = (double*) malloc (N * N * sizeof(double));
    double *phi = (double*) malloc (N * N * sizeof(double));
    double *a_x = (double*) malloc (N * N * sizeof(double));
    double *a_y = (double*) malloc (N * N * sizeof(double));
    double *particle_pos  = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_vel  = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_mass = (double*) malloc (N_p * sizeof(double));
    fftw_complex *rho_k   = (fftw_complex*) fftw_malloc (N * N * sizeof(fftw_complex));

    if (debug) {
        // random_particle_initialization(N_p, particle_pos, particle_vel, particle_mass, L);
        debug_initialization(N_p, particle_pos, particle_vel, particle_mass, L);
    } else {
        solar_sys_initialization(N_p, particle_pos, particle_vel, particle_mass, L);
    }

#ifdef MARSHAL
    Marshaller marshaller("particles.txt", L, N, N_p, particle_mass);
    marshaller.marshal(particle_pos);
#endif

    fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
    fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

    // TIME STEP ////////////////////////////////
    auto t_start = chrono::high_resolution_clock::now();

    for (int t=1; t<T; t++) {

        compute_rho(N_p, N, rho, particle_mass, particle_pos, delta_d);

        fftw_execute(rho_plan);

        compute_phi_k(N, L, rho_k, G);

        fftw_execute(phi_plan);

        // TODO: Scaling delta_t

        compute_accelerations(N, a_x, a_y, phi, delta_d);

        update_particles(N_p, N, particle_pos, particle_vel, a_x, a_y, delta_t, delta_d, L);

#ifdef MARSHAL
        marshaller.marshal(particle_pos);
#endif
    }

    auto t_end = chrono::high_resolution_clock::now();
    double total_time_ns = chrono::duration_cast<chrono::nanoseconds>(t_end-t_start).count();
    double average_time_ms = total_time_ns / (T * 1e6);
    cout << "Average time per step in milliseconds: " << average_time_ms << endl;

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
