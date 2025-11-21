// from protons on silicons, plot energy in silicon vs beam-like stopping point to find events with Ex=0 (Mauss thesis,
// figs 4.12-4.13)

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "./theoretical_pESil_vs_BPS.cxx"

void pESil_vs_BSP()
{
    // Read data
    auto fIn {"../Outputs/tree_pid_17F_p_p.root"};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {"PID_Tree", fIn};

    // Filter one silicon layer events
    auto df {d.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})};

    // Find beam-like stopping point (=RPx + TLheavy)
    df = df.Define("BSP",
                   [&](const ActRoot::MergerData& m)
                   {
                       double bsp = m.fRP.X() + m.fHeavy.fTL;
                       return bsp;
                   },
                   {"MergerData"});


    int thetamin {5};
    int thetamax {25};
    int step {5};
    std::vector<ROOT::TThreadedObject<TH2D>*> hsEpBSP;

    std::map<std::string, std::vector<TGraph*>> hsTheo;
    std::map<std::string, std::tuple<double, EColor, int, double, double>> theoArgs {
        {"^{17}F gs", {0.0, kRed, 1, 183., 12.3}},
        {"^{17}F* 0.495 MeV", {0.495, kBlack, 9, 166., 9.4}},
        {"^{17}F* 3.104 MeV", {3.104, kMagenta, 10, 186., 1.5}}};

    int idx {0};
    for(double theta = thetamin; theta < thetamax; theta += step)
    {
        hsEpBSP.push_back(new ROOT::TThreadedObject<TH2D>(
            TString::Format("hEpBSP%d", idx),
            TString::Format("#theta_{lab} = %.0f^{o}; Beam-like Stopping Point [cm]; E_{Si} [MeV]", theta), 65, 150,
            280, 150, 0, 15));

        for(const auto& [key, args] : theoArgs)
        {
            auto [Eex, color, mstyle, xlabel, ylabel] = args;
            hsTheo[key].push_back(calcTheo_pESil_vs_BSP(theta, Eex, color, mstyle));
        }

        idx++;
    }

    // Initialize slot 0 to not crash
    for(auto& h : hsEpBSP)
        h->GetAtSlot(0);

    // Fill histograms
    df.ForeachSlot(
        [&](unsigned int slot, double bsp, const ActRoot::MergerData& m)
        {
            // get the hstogram we have to fill
            for(size_t i = 0; i < hsEpBSP.size(); i++)
            {
                double theta = thetamin + i * step;
                if(m.fThetaLight >= theta && m.fThetaLight < theta + step)
                {
                    hsEpBSP[i]->GetAtSlot(slot)->Fill(bsp, m.fLight.fEs.front());
                }
            }
        },
        {"BSP", "MergerData"});

    // Draw
    auto* c0 {new TCanvas("c0", "ESil vs BSP", 900, 600)};
    c0->DivideSquare(hsEpBSP.size());
    int p {1};


    for(auto& h : hsEpBSP)
    {
        c0->cd(p);
        h->Merge()->DrawClone("colz");

        // Draw theretical graphs
        for(const auto& [key, g] : hsTheo)
        {
            g[p - 1]->DrawClone("same");
        }
        // Draw labels
        for(const auto& [key, args] : theoArgs)
        {
            auto [Eex, color, mstyle, xlabel, ylabel] = args;
            auto t = new TLatex(xlabel, ylabel,key.c_str());
            t->SetTextColor(color);
            t->Draw("same");
        }

        p++;
    }
}