#pragma once

namespace S2029
{

// Global parameters used throughout the S2029 analysis

// Pressure calculated comparing to drift velocity alphas in S2029/Calibrations/Actar/fineTuneSRIMfromAlphas.cxx
inline constexpr double pressure {765.};//{775.}; // mbar

// Initial Beam Energy calculated in S2029/Macros/Beam/calcEBeamIni.cxx
inline constexpr double EBeamIni {3.87};//{3.905}; // MeV/u

// Total beam particles on target from S2029/PostAnalysis/Pipes/Pipe_Beam.cxx
inline constexpr double CFADiv {176922.}; 

// Target density from LISE++
inline constexpr double Nt {1.25e21}; // atoms/cm2 for 256mm target thickness at 775mbar

} // namespace S2029