// pESil_vs_BSP events have a few poorly reconstructed events on the left side of the band. I want to see if they
// correspond to a specific silicon that may have incorrect calibration.

#include "ActCutsManager.h"
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

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

void investigate_pESil_vs_BSP()
{
    // Read data
    auto fIn {"./Outputs/tree_pESil_vs_BSP_good_events.root"};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {"PID_Tree", fIn};


    // read cut on poorly reconstructed events
    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("cut", TString::Format("./Outputs/pESil_vs_BSP_cut_good_events.root").Data());

    auto gated {d.Filter(
        [&](double bsp, ActRoot::MergerData& m)
        {
            if(cuts.GetCut("cut"))
                return cuts.IsInside("cut", bsp, m.fLight.fEs.front());
            else
                return false;
        },
        {"BSP", "MergerData"})};

    std::array<int, 12> f0 {}, l0 {}, r0 {};


    gated.Foreach(
        [&](ActRoot::MergerData& m)
        {
            auto layer {m.fLight.fLayers.front()};
            auto n {m.fLight.fNs.front()};
            if(layer == "f0")
                f0[n] += 1;
            if(layer == "r0")
                r0[n] += 1;
            if(layer == "l0")
                l0[n] += 1;
        },
        {"MergerData"});


    for(const auto& layer : {std::pair {"f0", f0}, std::pair {"l0", l0}, std::pair {"r0", r0}})
    {
        std::cout << layer.first << std::endl;
        for(std::size_t i = 0; i < layer.second.size(); i++)
            std::cout << i << " " << layer.second[i] << " " << std::endl;
    }
}