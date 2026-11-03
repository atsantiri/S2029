#pragma once

namespace S2029
{

// Global parameters used throughout the S2029 analysis

// Pressure calculated comparing to drift velocity alphas in S2029/Calibrations/Actar/fineTuneSRIMfromAlphas.cxx
inline constexpr double pressure {765.}; // mbar

// Initial Beam Energy calculated in S2029/Macros/Beam/calcEBeamIni.cxx
inline constexpr double EBeamIni {3.87}; // MeV/u
// For 770 mbar-> EbeamIni = 3.89+-0.03 MeV/u
// For 765 mbar-> EbeamIni = 3.87+-0.03 MeV/u
// For 760 mbar-> EbeamIni = 3.85+-0.03 MeV
// For 755 mbar-> EbeamIni = 3.83+-0.03 MeV
// For 740 mbar-> EbeamIni = 3.78+-0.03 MeV
// For 730 mbar-> EbeamIni = 3.75+-0.03 MeV


// Total beam particles on target from S2029/PostAnalysis/Pipes/Pipe_Beam.cxx
inline constexpr double CFADiv {176922.}; 

// Target density from LISE++
inline constexpr double Nt {1.25e21}; // atoms/cm2 for 256mm target thickness at 775mbar

} // namespace S2029