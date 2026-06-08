#include "ROOT/RDataFrame.hxx"
#include "AngFitter.h"
#include "AngGlobals.h"
#include "AngIntervals.h"
#include "FitInterface.h"
#include "Interpolators.h"
#include "../HistConfig.h"

void normalizeEx()
{
    // Pipe 3 output
    ROOT::RDataFrame df{"Final_Tree", "../Outputs/tree_ex_17F_p_p_3.84.root"};

    // Simulation output
    ROOT::RDataFrame dsim{"SimulationTTree", "../Input/Simu_1H_760mbar.root"};

    // Book histograms
    auto hEx{df.Histo1D({"hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (4. - (-2.)) / 200 * 1e3), 200, -2, 4}, "RecEx")};
    auto hCM{df.Histo2D({"hCM", "CM;#theta_{CM};E [MeV]", 300, 0, 120, 300, 0, 60}, "RecThetaCM", "EVertex")};

    // Init intervals
    double thetaMin{40};
    double thetaMax{80};
    double thetaStep{10};
    Angular::Intervals ivs{thetaMin, thetaMax, HistConfig::Ex, thetaStep, 0};
    df.Foreach([&](double thetacm, double ex)
               { ivs.Fill(thetacm, ex); }, {"RecThetaCM", "RecEx"});
    ivs.Draw();

    // Efficiency
    Interpolators::Efficiency eff;
    eff.Add("g0", "../Input/Simu_1H_760mbar.root", "effCM");
    // Draw to check is fine
    eff.Draw();
}