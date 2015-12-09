#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cmath>
#include <complex.h>
#include <fftw3.h>

using namespace std;

int main() {
    int N=128, N_p = 3;
    double L = 100.0;
    double *rho = (double*) malloc (N * N * sizeof(double));
    double *E_x = (double*) malloc (N * N * sizeof(double));
    double *E_y = (double*) malloc (N * N * sizeof(double));
    fftw_complex *rho_k_x = (fftw_complex*) fftw_malloc (N * N * sizeof(fftw_complex));
    fftw_complex *rho_k_y = (fftw_complex*) fftw_malloc (N * N * sizeof(fftw_complex));
    double *particle_pos = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_vel = (double*) malloc (N_p * 2 * sizeof(double));
    double *particle_mass = (double*) malloc (N_p * sizeof(double));
    double delta_t = 0.1;
    double delta_d = L/N;


    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> pos_dis(0.0, L);
    uniform_real_distribution<> mass_dis(0.1, 3.0);
    ofstream fp;
    fp.open("initial.csv");
    for (int i=0; i<N_p; i++) {
        particle_pos[2*i]   = pos_dis(gen);
        particle_pos[2*i+1] = pos_dis(gen);
        fp << particle_pos[2*i] << " " << particle_pos[2*i+1] << endl;
        particle_vel[2*i]   = 0.0;
        particle_vel[2*i+1] = 0.0;
        particle_mass[i]    = mass_dis(gen);
    }
    fp.close();
    
    fftw_plan rho_plan =  fftw_plan_dft_r2c_2d(N, N, rho, rho_k_x, FFTW_ESTIMATE);
    fftw_plan E_x_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k_x, E_x, FFTW_ESTIMATE);
    fftw_plan E_y_plan =  fftw_plan_dft_c2r_2d(N, N, rho_k_y, E_y, FFTW_ESTIMATE);

    double time = 0.0;
    for (int t=1; t<30; t++) {

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
            for (int i=0; i<N; i++) {
                if (i==0 && j==0) continue;
                

                if (i >= N/2)  n_i = (i - N); // Negative frequency mapping in x
                else n_i = i;

                if (j >= N/2)  n_j = (j - N); // Negative frequency mapping in y
                else n_j = j;

                kx_i = n_i * 2 * M_PI / L; // Spatial frequency in x
                ky_j = n_j * 2 * M_PI / L; // Spatial frequency in y
                k_sq = kx_i * kx_i + ky_j * ky_j;
                mult = -4 * M_PI / k_sq;
                rho_k_y[j*N+i][0] = rho_k_x[j*N+i][0];
                rho_k_y[j*N+i][1] = rho_k_x[j*N+i][1];

                a = rho_k_x[j*N+i][0] * mult * kx_i;
                rho_k_x[j*N+i][0] = - rho_k_x[j*N+i][1] * mult;
                rho_k_x[j*N+i][1] = a;

                a = rho_k_y[j*N+i][0] * mult * ky_j;
                rho_k_y[j*N+i][0] = - rho_k_y[j*N+i][1] * mult;
                rho_k_y[j*N+i][1] = a;
            }
        }
        fftw_execute(E_x_plan);
        fftw_execute(E_y_plan);

        char fname[10];
        sprintf(fname, "dt%d.csv", t);

        fp.open(fname);
        for (int i=0; i<N_p; i++) {
            ind_x =  floor(particle_pos[2*i]/delta_d);
            ind_y =  floor(particle_pos[2*i+1]/delta_d);
            if (ind_x >= N || ind_x < 0 || ind_y >= N || ind_y < 0) {
                fp << particle_pos[2*i] << " " << particle_pos[2*i+1] << endl;
                continue;
            }
            particle_vel[2*i] += E_x[ind_y*N + ind_x]/ (N*N) * delta_t;
            particle_vel[2*i+1] += E_y[ind_y*N + ind_x]/ (N*N) * delta_t;
            particle_pos[2*i] += particle_vel[2*i] * delta_t;
            particle_pos[2*i+1] += particle_vel[2*i+1] * delta_t;
            fp << particle_pos[2*i] << " " << particle_pos[2*i+1] << endl;
        }
        fp.close();
    }


    fftw_destroy_plan(E_x_plan);
    fftw_destroy_plan(E_y_plan);
    fftw_destroy_plan(rho_plan);
    free(particle_pos);
    free(particle_vel);
    free(particle_mass);
    free(rho);
    free(E_x);
    free(E_y);
    fftw_free(rho_k_x);
    fftw_free(rho_k_y);
}
