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
    "  - o -- number of OMP threads (24)\n"
    "  - s -- number of time steps (1000)\n";

void strong_scaling(int t_steps, double* particle_pos, double* particle_vel,
                    double* particle_mass, int N_p, int N,
                    double L, double* a_x, double* a_y, double* rho, double* phi,
                    fftw_complex* rho_k, double delta_d, double delta_t) {
    FILE *fp;
    fp = fopen("strong_scaling.csv", "w+");
    double t_start, t_end, t_threadrun;
    int thread_max = 26, i;

    for (i = 1; i <= thread_max; i++) {
        omp_set_num_threads(i);
        random_particle_initialization(N_p, L,
                                       particle_vel, particle_pos,
                                       particle_mass);

        fftw_plan_with_nthreads(i);
        fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
        fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

        t_threadrun = 0.0;

        for (int t=1; t<t_steps; t++) {
            t_start = omp_get_wtime();

            compute_rho(N_p, N, delta_d, rho,
                        particle_mass, particle_pos);

            fftw_execute(rho_plan);

            compute_phi_k(N, L, rho_k);

            fftw_execute(phi_plan);

            // TODO: Scaling delta_t

            compute_accelerations(N, a_x, a_y, phi, delta_d);

            update_particles(N_p, N, delta_t, delta_d, L,
                             particle_pos, particle_vel, a_x, a_y);


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
    int num_threads = 24;

    int error_code = fftw_init_threads();
    if (error_code == 0) {
        cout << "Error initializing FFTW threads. Exiting..." << endl;
        exit( EXIT_FAILURE );
    }

    // Option processing
    extern char* optarg;
    const char* optstring = "hp:n:l:t:s:o:";
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
        case 'o': num_threads = atoi(optarg);  break;
        }
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

    random_particle_initialization(N_p, L,
                                   particle_vel, particle_pos,
                                   particle_mass);

#ifdef MARSHAL
    Marshaller marshaller("particles.txt", L, N, N_p, particle_mass);
    marshaller.marshal(particle_pos);
#endif

    omp_set_num_threads(num_threads);
    fftw_plan_with_nthreads(num_threads);
    fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k, FFTW_MEASURE);
    fftw_plan phi_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k, phi, FFTW_MEASURE);

    // TIME STEP ////////////////////////////////
    double t_start = omp_get_wtime();

    #pragma omp parallel
    for (int t=1; t<t_steps; t++) {
        // #pragma omp single
        // compute_rho(N_p, N, delta_d, rho,
        //             particle_mass, particle_pos);

        {
            // Zero-out rho for reuse
            #pragma omp for
            for (int i = 0; i < N * N; ++i) {
                rho[i] = 0.0;
            }

            int ind_x, ind_y, index, mass;
            #pragma omp for
            for (int i=0; i<N_p; i++) {
                ind_x =  floor(particle_pos[2*i]/delta_d);
                ind_y =  floor(particle_pos[2*i+1]/delta_d);
                index = ind_y*N + ind_x;
                mass = particle_mass[i];
                #pragma omp critical
                rho[index] += mass;
            }

            // Scale to represent density
            double scaling_factor = delta_d * delta_d;
            #pragma omp for
            for (int i=0; i<N*N; i++) {
                rho[i] /= scaling_factor;
            }
        }

        #pragma omp single
        fftw_execute(rho_plan);

        // #pragma omp single
        // compute_phi_k(N, L, rho_k);
        {
            double n_i, n_j, kx_i, ky_j, k_sq, scaling_factor;

            #pragma omp for
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

        #pragma omp single
        fftw_execute(phi_plan);

        // TODO: Scaling delta_t

        // #pragma omp single
        // compute_accelerations(N, a_x, a_y, phi, delta_d);
        {
            // 2 delta_d for finite difference and N*N for FFTW normalization
            double scaling_factor = 2 * delta_d * N * N;
            #pragma omp for collapse(2) nowait
            for (int j=0; j<N; j++) {
                for (int i=1; i<N-1; i++) {
                    a_x[j*N + i] = (-phi[j*N + i-1] +  phi[j*N + i+1]) / scaling_factor;
                }
            }

            #pragma omp for collapse(2) nowait
            for (int j=1; j<N-1; j++) {
                for (int i=0; i<N; i++) {
                    a_y[j*N + i] = (-phi[(j-1)*N + i] +  phi[(j+1)*N + i]) / scaling_factor;
                }
            }

            // Boundary conditions for phi are periodic
            // Left-right edges
            #pragma omp for nowait
            for (int j=0; j<N; j++) {
                a_x[j*N] = (-phi[j*N + (N-1)] + phi[j*N + 1]) / scaling_factor;
                a_x[j*N+(N-1)] = (-phi[j*N + N-2] + phi[j*N]) / scaling_factor;
            }

            // Top-bottom edges
            #pragma omp for
            for (int i=0; i<N; i++) {
                a_y[i] = (-phi[(N-1)*N +i] +  phi[N + i]) / scaling_factor;
                a_y[(N-1)*N + i] = (-phi[(N-2)*N + i] + phi[i])/ scaling_factor;
            }
        }

        // #pragma omp single
        // update_particles(N_p, N, delta_t, delta_d, L,
                         // particle_pos, particle_vel, a_x, a_y);
        {
            int ind_x, ind_y;
            #pragma omp for
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

#ifdef MARSHAL
        #pragma omp single
        marshaller.marshal(particle_pos);
#endif
    }
    double t_end = omp_get_wtime();
    double total_time_ms = t_end - t_start;
    double average_time_ms = total_time_ms / t_steps;

    cout << "ppm_omp_full "
         << L << " "
         << N << " "
         << N_p << " "
         << delta_t << " "
         << t_steps << " "
         << num_threads << " "
         << total_time_ms << " "
         << average_time_ms << " "
         << endl;

    // strong_scaling(t_steps, particle_pos, particle_vel, particle_mass, N_p, N,
            // L, a_x, a_y, rho, phi, rho_k, delta_d, delta_t);

    fftw_destroy_plan(rho_plan);
    fftw_destroy_plan(phi_plan);
    free(particle_pos);
    free(particle_vel);
    free(particle_mass);
    free(rho);
    free(a_x);
    free(a_y);
    fftw_free(rho_k);
    fftw_cleanup_threads();
}
