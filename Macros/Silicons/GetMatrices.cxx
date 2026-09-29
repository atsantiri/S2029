#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"

#include <ios>
#include <iostream>
#include <vector>

#include "../../PostAnalysis/Utils.cxx"

void GetMatrices(TString mode)
{
    mode.ToLower();
    std::cout << "Running for : " << mode << '\n';

    // Assign things depending on mode
    std::vector<int> idxs {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    bool isSide {mode.Contains("f") ? false : true};
    std::cout << "isSide ? " << std::boolalpha << isSide << '\n';

    // Read data
    ROOT::EnableImplicitMT();
    ActRoot::DataManager datman {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {datman.GetChain()};
    ROOT::RDataFrame df {*chain};

    // Filter
    auto gated {df.Filter(
        [mode](ActRoot::MergerData& m)
        {
            if(mode.Contains("l"))
                return S2029::isLeft(m);
            if(mode.Contains("r"))
                return S2029::isRight(m);
            if(mode.Contains("f"))
                return S2029::isFront(m);
            return false;
        },
        {"MergerData"})};

    // Book histograms
    int xybins {300};
    std::pair<double, double> xlims {-20, 430};
    int zbins {300};
    std::pair<double, double> zlims {-20, 430};
    auto hSP {gated.Histo2D(
        {"hSP", "SP;X | Y [mm];Z [mm]", xybins, xlims.first, xlims.second, zbins, zlims.first, zlims.second},
        isSide ? "fSP.fCoordinates.fX" : "fSP.fCoordinates.fY", "fSP.fCoordinates.fZ")};

    // Gate and get projections
    std::map<int, ROOT::TThreadedObject<TH1D>> pxys, pzs;
    for(const auto& idx : idxs)
    {
        pxys.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                     std::forward_as_tuple(TString::Format("pxy%d", idx),
                                           TString::Format("X or Y proj %d;X | Y [mm]", idx), xybins, xlims.first,
                                           xlims.second));
        pzs.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                    std::forward_as_tuple(TString::Format("pz%d", idx), TString::Format("Z proj %d;Z [mm]", idx), zbins,
                                          zlims.first, zlims.second));
    }
    gated.ForeachSlot(
        [&](unsigned int slot, ActRoot::MergerData& data)
        {
            auto idx {data.fSilNs.front()};
            if(pxys.count(idx))
            {
                pxys[idx].GetAtSlot(slot)->Fill(isSide ? data.fSP.X() : data.fSP.Y());
                pzs[idx].GetAtSlot(slot)->Fill(data.fSP.Z());
            }
        },
        {"MergerData"});
    // Merge
    for(auto map : {&pxys, &pzs})
        for(auto& [_, h] : *map)
            h.Merge();

    // plot
    auto* c1 {new TCanvas("c1", ("Hists canvas for " + mode).Data())};
    hSP->DrawClone("colz");

    auto* cpxy {new TCanvas("cpxy", "X | Y projection canvas")};
    cpxy->DivideSquare(pxys.size());
    int pad {0};
    for(auto& [_, p] : pxys)
    {
        cpxy->cd(pad + 1);
        p.GetAtSlot(0)->DrawClone();
        pad++;
    }
    auto* cpz {new TCanvas("cpz", "Z projection canvas")};
    cpz->DivideSquare(pzs.size());
    pad = 0;
    for(auto& [_, p] : pzs)
    {
        cpz->cd(pad + 1);
        p.GetAtSlot(0)->DrawClone();
        pad++;
    }
    // Write them
    auto fout {std::make_unique<TFile>(TString::Format("../Outputs/%s_histograms.root", mode.Data()), "recreate")};
    fout->cd();
    for(auto& [_, h] : pxys)
        h.GetAtSlot(0)->Write();
    for(auto& [_, h] : pzs)
        h.GetAtSlot(0)->Write();
}
