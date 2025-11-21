#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"
#include <ROOT/TThreadedObject.hxx>


#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe2_Ex(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {"PID_Tree", filename};

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
    std::string srimName {};
    if(light == "p")
        srimName = "1H";
    else if(light == "4He")
        srimName = "4He";
    srim->ReadTable(light,
                    TString::Format("../Simulation/SRIM/%s_H2-iC4H10_95-5_760mbar.txt", srimName.c_str()).Data());
    srim->ReadTable(beam, TString::Format("../Simulation/SRIM/%s_H2-iC4H10_95-5_760mbar.txt", beam.c_str()).Data());
    srim->ReadTable("heavy", "../Simulation/SRIM/14O_H2-iC4H10_95-5_760mbar.txt");

    // Init particles
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};

    // Initial energy of beam at pad plane entrance
    double EBeamIni {3.84}; // MeV/u


    // Find beam-like stopping point (=RPx + TLheavy)
    auto df {d.Define("EBeam",
                      [&](const ActRoot::MergerData& d)
                      {
                          auto ret {srim->Slow(beam, EBeamIni * pb.GetAMU(),
                                               d.fRP.X())}; // here RPx comes from merger data so it's in mm, not pads
                          if(ret <= 0)
                              ret = 1111111;
                          return ret;
                      },
                      {"MergerData"})
                 .Define("BSP",
                         [&](const ActRoot::MergerData& m)
                         {
                             double bsp = m.fRP.X() + m.fHeavy.fTL;
                             return bsp;
                         },
                         {"MergerData"})};

    int thetamin {5};
    int thetamax {40};
    int step {5};
    std::vector<ROOT::TThreadedObject<TH2D>*> hsEpBSP;

    int idx {0};
    for(double theta = thetamin; theta < thetamax; theta += step)
    {
        hsEpBSP.push_back(new ROOT::TThreadedObject<TH2D>(
            TString::Format("hEpBSP%d", idx),
            TString::Format("#theta_{lab} [%.2f, %.2f); Beam-like Stopping Point [cm]; E_{Si} [MeV]", theta, theta + step), 65, 150, 280, 150,
            0, 15));
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
    auto* c0 {new TCanvas("c0", "ESil vs BSP")};
    c0->DivideSquare(hsEpBSP.size());
    int p {1};
    for(auto& h : hsEpBSP)
    {
        c0->cd(p);
        h->Merge()->DrawClone("colz");
        p++;
    }
                         
}
#endif