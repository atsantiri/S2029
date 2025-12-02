#ifndef Pipe2_ESil_BSP_cxx
#define Pipe2_ESil_BSP_cxx

#include "ActCutsManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"

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

#include "../Macros/theoretical_pESil_vs_BPS.cxx"


void Pipe2_ESil_BSP(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {"PID_Tree", filename};

    // Init particles
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};

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
    int thetamax {35};
    int step {5};

    // make histos
    std::vector<ROOT::TThreadedObject<TH2D>*> hsEpBSP;
    std::map<std::string, std::vector<TGraph*>> hsTheo;
    std::map<std::string, std::tuple<double, EColor, int, double, double>> theoArgs {
        {"^{17}F gs", {0.0, kRed, 1, 183., 12.3}}, {"^{17}F* 0.495 MeV", {0.495, kBlack, 9, 180., 2.}}};
    // {"^{17}F* 3.104 MeV", {3.104, kMagenta, 10, 186., 1.5}}}; // {label, Eex, linecolor, linestyle, xlabel, ylabel}

    int idx {0};
    for(double theta = thetamin; theta < thetamax; theta += step)
    {
        hsEpBSP.push_back(new ROOT::TThreadedObject<TH2D>(
            TString::Format("hEpBSP%d", idx),
            TString::Format("#theta_{lab} [%.0f^{o}, %.0f^{o}); Beam-like Stopping Point [cm]; E_{Si} [MeV]", theta,
                            theta + step),
            65, 150, 280, 150, 0, 15));

        for(const auto& [key, args] : theoArgs)
        {
            auto [Eex, color, mstyle, xlabel, ylabel] = args;
            hsTheo[key].push_back(calcTheo_pESil_vs_BSP(theta + step / 2, Eex, color, mstyle));
        }

        idx++;
    }

    // Initialize slot 0 to not crash
    for(auto& h : hsEpBSP)
        h->GetAtSlot(0);

    // Fill exp histograms
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

    int i = 0;
    for(auto& h : hsEpBSP)
    {
        c0->cd(i + 1);
        h->Merge()->DrawClone("colz");

        for(const auto& [key, graphs] : hsTheo)
            graphs[i]->DrawClone("same");

        // Draw labels
        for(const auto& [key, args] : theoArgs)
        {
            auto [Eex, color, mstyle, xlabel, ylabel] = args;
            auto t = new TLatex(xlabel, ylabel, key.c_str());
            t->SetTextColor(color);
            t->Draw("same");
        }
        i++;
    }


    // If cut on Ex=0 band present, apply
    ActRoot::CutsManager<std::string> cuts;
    // cuts.ReadCut("cut", TString::Format("./Cuts/ESil_BSP_%s_%s.root", light.c_str(), beam.c_str()).Data());

    auto listOfCuts {cuts.GetListOfKeys()};
    if(listOfCuts.size())
    {
        // Apply cut and save in file
        auto gated {df.Filter(
            [&](double bsp, ActRoot::MergerData& m)
            {
                if(cuts.GetCut("cut"))
                    return cuts.IsInside("cut", bsp, m.fLight.fEs.front());
                else
                    return false;
            },
            {"BSP", "MergerData"})};
        auto name {
            TString::Format("./Outputs/tree_ESil_BSP_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
        std::cout << "Saving PID_Tree with Ex=0 in file : " << name << '\n';
        gated.Snapshot("PID_Tree", name.Data());

        // Draw the cut
        for(int i = 0; i < hsEpBSP.size(); i++)
        {
            c0->cd(i + 1);
            cuts.DrawCut("cut");
        }
    }
    else
    { // there is no cut but i still want the file
        auto name {
            TString::Format("./Outputs/tree_ESil_BSP_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
        std::cout << "Saving PID_Tree with Ex=0 in file : " << name << '\n';
        df.Snapshot("PID_Tree", name.Data());
    }
}
#endif