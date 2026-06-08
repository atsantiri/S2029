#include "TString.h"

#include <stdexcept>
#include <string>
#include <vector>

// #include "./Plotter.cxx"
#include "./Simulation_S2029.cxx"

// what
// if simu = runs simulation
// if plot = plots results
// standalone (only applies to simu setting):
// if true, runs only first item in Ex vector and plots in-simulation results
// if false, runs all Ex simulations but doesn't plot
void Runner(TString what = "simu", bool standalone = true)
{
  // Settings
  // Names of particles
  std::string beam{"17F"};
  std::string target{"1H"};
  std::string light{"1H"};
  std::string heavy{"17F"};
  // Phase space reactions: when the heavy decays by proton or neutron emission
  // So we have something like: 4He + n + 17N (needs to be simulated to be
  // included as background in fits)
  // int neutronPS{0};  // number of neutrons in final state
  // int protonPS{0};   // number of protons in final state
  double T1{3.84};   // Beam energy at entrance of pad plane
  int pressure{760}; // mbar

  if (what.Contains("simu"))
    Simulation_S2029(beam, target, light, heavy, T1, 0, pressure, standalone);
  if (what.Contains("plot"))
  {
    throw std::runtime_error("I'll work on the plotter, maybe");
    // Plotter({0.}, beam, target, light, T1, pressure);
  }
}
