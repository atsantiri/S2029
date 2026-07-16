#include "ROOT/RDataFrame.hxx"
#include "AngFitter.h"
#include "AngGlobals.h"
#include "AngIntervals.h"
#include "FitInterface.h"
#include "Interpolators.h"
#include "../HistConfig.h"

#include "TLine.h"

void normalizeEx()
{
    // Pipe 3 output
    ROOT::RDataFrame df{"Final_Tree", "../Outputs/tree_ex_17F_p_p_3.84_sil.root"};

    // Simulation output
    TString simuFile = "../Input/Simu_1H_760mbar.root";
    ROOT::RDataFrame dsim{"SimulationTTree", simuFile};

    // Book histograms
    auto hEx{df.Histo1D({"hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (4. - (-2.)) / 200 * 1e3), 200, -2, 4}, "RecEx")};
    auto hCM{df.Histo2D(HistConfig::EcnThetaCM, "RecThetaCM", "RecECN")};

    // Init intervals
    double thetaMin{40};
    double thetaMax{170};
    double thetaStep{10};
    Angular::Intervals ivsEx{thetaMin, thetaMax, HistConfig::Ex, thetaStep, 0};
    df.Foreach([&](double thetacm, double ex)
               { ivsEx.Fill(thetacm, ex); }, {"RecThetaCM", "RecEx"});
    ivsEx.Draw();

    Angular::Intervals ivsEcm{thetaMin, thetaMax, HistConfig::ECN, thetaStep, 0};
    df.Foreach([&](double thetacm, double ecn)
               { ivsEcm.Fill(thetacm, ecn); }, {"RecThetaCM", "RecECN"});
    ivsEcm.Draw();
    TCanvas *c1 = gPad->GetCanvas();
    // Overlay known states
    std::vector<double> states{{4519, 4561, 4590, 5090, 5146, 5453, 6150, 6297, 6353, 7059, 7350, 7713, 7910, 7950}};
    double qval{3.923}; // qvalue of 17F+p

    for (auto &s : states)
    {
        TLine *st = new TLine(s*1e-3, 0, s*1e-3, 150);
        st->SetLineColorAlpha(16, 0.2);
        st->SetLineStyle(2);
        for (int i = 1; i <= ivsEcm.GetHistos().size(); i++)
        {
            c1->cd(i);
            st->Draw("same");
        }
    }

    // Init fitter
    // Angular::Fitter fitter{&ivs};
    // fitter.SetAllowFreeMean(false);
    // // fitter.SetFreeMeanRange(0.1);
    // fitter.Configure("./Outputs/fit_pp.root");
    // fitter.Run();
    // fitter.Draw();
    // fitter.DrawCounts();

    // Efficiency
    Interpolators::Efficiency eff;
    eff.Add("g0", simuFile.Data(), "effCM");
    // Draw to check is fine
    eff.Draw();


}