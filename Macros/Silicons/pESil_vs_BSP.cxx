// from protons on silicons, plot energy in silicon vs beam-like stopping point to find events with Ex=0 (Mauss thesis, figs 4.12-4.13)
//
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


void pESil_vs_BSP()
{
    // Read data
    auto fIn {"../../PostAnalysis/Outputs/tree_pid_17F_p_p.root"};
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


    int thetamin {20};
    int thetamax {100};
    int step {20};
    std::vector<ROOT::TThreadedObject<TH2D>*> hsEpBSP;

    int idx {0};
    for(double theta = thetamin; theta < thetamax; theta += step)
    {
        hsEpBSP.push_back(new ROOT::TThreadedObject<TH2D>(
            TString::Format("hEpBSP%d", idx),
            TString::Format("#theta_{lab} [%.2f, %.2f); BSP [cm]; E_{Si} [MeV]", theta, theta + step), 150, 0, 256, 150,
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