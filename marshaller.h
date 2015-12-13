#ifndef MARSHALLER_H
#define MARSHALLER_H

#include <fstream>

// During each step of a particle simulation, we update the position and
// velocity of each particle. After each step, we write the position of each
// particle to a file: a process known as marshalling. These text files can
// then be parsed by other programs to generates visualizations of the
// simulation. The Marshal class helps marshal simulation data to a file in a
// convenient and consistent way.
class Marshaller {
  public:
    // Construct a Marshaller for an L by L simulation with N by N grid points
    // and N_p particles and their corresponding masses (mass). All data is
    // written to `filename`.
    Marshaller(const std::string& filename, const double L, const int N,
               const int N_p, const double *mass);

    Marshaller(const Marshaller&) = delete;

    ~Marshaller();

    // `marhsal` receives two arrays of size `N_p`. The second list is a list
    // of particle positions; some of the particles are within the L by L
    // simulation space while some are outside of it. The first list specifies
    // which positions are within the space.
    void marshal(const bool *valid, const double *positions);

  private:
    std::ofstream f_; // underlying file
    const int N_p_;   // # of particles
    int n_;           // # of time steps == # of times marshal is called
};

#endif // MARSHALLER_H
