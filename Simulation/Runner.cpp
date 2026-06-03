#include "TString.h"

#include <stdexcept>
#include <string>
#include <vector>

// #include "./Plotter.cpp"
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
  int neutronPS{0};  // number of neutrons in final state
  int protonPS{0};   // number of protons in final state
  double T1{};       // Beam energy at entrance of pad plane // 1.525 with 1.8 //1.486
                     // with 1.774
  int pressure{760}; // mbar
  T1 = 3.84;         // MeV/u

  std::vector<double> Eexs;
  if (neutronPS == 0 && protonPS == 0)
    Eexs = {0};
  else if (neutronPS > 0 && protonPS == 0)
    Eexs = {0}; // only gs for n phase space
  else if (neutronPS == 0 && protonPS > 0)
    Eexs = {0};
  else
    throw std::runtime_error(
        "No confs with neutronPS and protonPS enabled at the same time");

  if (what.Contains("simu"))
  {
    for (const auto &Eex : Eexs)
    {
      Simulation_S2029(beam, target, light, heavy, neutronPS, protonPS, T1, Eex,
                       pressure, standalone);
      if (standalone)
        break;
    }
  }
  if (what.Contains("plot"))
  {
    throw std::runtime_error(
        "Runner for E837 hasn't implemented yet a Plotter function");
    // Plotter(Eexs, beam, target, light, heavy, T1, neutronPS, protonPS);
  }
}
