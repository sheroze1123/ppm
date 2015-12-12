#ifndef MARSHALLER_H
#define MARSHALLER_H

#include <fstream>

// During each step of a particle simulation, we update the position and
// velocity of each particle. After each step, we write the position and
// velocity of each particle to a file, and another program parses the file and
// generates a visualization of the simulation. The Marshal class helps write
// the values to a file in a convenient and consistent way.
class Marshaller {
  public:
    // Construct a Marshaller for an L by L simulation with N by N grid points
    // and N_p particles and their corresponding masses (mass). All data is
    // written to `filename`.
    Marshaller(const std::string& filename, const double L, const int N,
               const int N_p, const double *mass);

    Marshaller(const Marshaller&) = delete;

    ~Marshaller();

    // Write the position and velocity of all N_p particles to the file managed
    // by this Marshaller.
    void marshal(const bool *valid, const double *positions);

  private:
    std::ofstream f_;
    const int N_p_; // number of particles
    int n_;         // number of times marshal is called
};

#endif // MARSHALLER_H
