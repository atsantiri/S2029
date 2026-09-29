#pragma once

namespace S2029
{

// Global variables used throughout the S2029 analysis

// Pressure calculated comparing to drift velocity alphas in S2029/Calibrations/Actar/fineTuneSRIMfromAlphas.cxx
inline constexpr double pressure {775.}; // mbar

// Initial Beam Energy calculated in S2029/Macros/Beam/calcEBeamIni.cxx
inline constexpr double EBeamIni {3.905}; // MeV/u

} // namespace S2029