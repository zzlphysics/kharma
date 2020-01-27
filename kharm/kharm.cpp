/*
 * K/HARM -- Implementation of the HARM scheme for GRMHD,
 * in C++ with Kokkos performance portability library
 *
 * Ben Prather
 */

#include "decs.hpp"
#include "diffuse.hpp"
#include "self_init.hpp"
#include "grid.hpp"

#if USE_MPI
// This is very not supported right now.
// Question whether we even use Boost MPI or stick to native bindings
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;
#endif

#include <cmath>
#include <iostream>
#include <random>
#include <sstream>

using namespace Kokkos;
using namespace std;

int main(int argc, char **argv)
{
  size_t ng = 3; // TODO add to dumps/take

#if USE_MPI
  mpi::environment env(argc, argv);
  mpi::communicator world;
#endif
  Kokkos::initialize(argc, argv);
  {
    std::cerr << "K/HARM v.alpha" << std::endl;
    std::cerr << "Using Kokkos environment:" << std::endl;
    DefaultExecutionSpace::print_configuration(std::cerr);
    std::cerr << std::endl;

    // TODO make right for parallel HDF5 later
    // HighFive::File input(argv[1], HighFive::File::ReadOnly);
    // auto prims_shape = input.getDataSet("/prims").getDimensions();
    // int n1 = prims_shape[0];
    // int n2 = prims_shape[1];
    // int n3 = prims_shape[2];
    // int nprim = prims_shape[3];
    // std::cout << "Input size: " << n1 << "x" << n2 << "x" << n3 << "x" << nprim << std::endl;

    // Contiguous array of just n1xn2xn3 for input
    // GridVarsHost h_prims_input("prims_input", n1, n2, n3, nprim);
    // input.getDataSet("/prims").read(h_prims_input.data());

    Grid G({128, 128, 128}, {0,0,0}, {1,1,1}, 3, 8);
    cerr << "Grid init" << std::endl;

    GridVarsHost h_vars_input = mhdmodes(G, 0);
    cerr << "Vars init" << std::endl;

    GridVars vars("all_vars", G.gn1, G.gn2, G.gn3, G.nprim);
    GridVars vars_temp("all_vars_temp", G.gn1, G.gn2, G.gn3, G.nprim);
    auto h_vars = create_mirror_view(vars);
    auto h_vars_temp = create_mirror_view(vars_temp);

    // Copy input (no ghosts, Host order) into working array (ghosts, device order)
    // deep_copy would do this automatically if not for ghosts (TODO try that?)
    parallel_for("diff_all", *(G.bulk_0),
                 KOKKOS_LAMBDA(int i, int j, int k) {
                   for (int p = 0; p < G.nprim; ++p)
                     h_vars(i + ng, j + ng, k + ng, p) = h_vars_input(i, j, k, p);
                 });

    // copy TO DEVICE
    deep_copy(vars, h_vars);

    for (int iter = 0; iter < 1000; iter++)
    {
      if (iter % 2 == 0)
        diffuse_all(vars, vars_temp);
      else
        diffuse_all(vars_temp, vars);
    }

    deep_copy(h_vars, vars);
  }
  Kokkos::finalize(); // Kokkos

  return 0;
}
