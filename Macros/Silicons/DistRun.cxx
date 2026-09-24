#include "ActDataManager.h"
#include "ActLine.h"
#include "ActMergerData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TFile.h"
#include "TH1.h"
#include "TStopwatch.h"
#include "TString.h"

#include "Math/Point3Dfwd.h"

#include <iostream>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "./Utils.cxx"

void DistRun(TString mode)
{
    mode.ToLower();
    bool isSide {mode.Contains("f") ? false : true};
    std::cout << "isSide ? " << std::boolalpha << isSide << '\n';

    ROOT::EnableImplicitMT();

    ActRoot::DataManager data {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {data.GetJoinedData()};

    ROOT::RDataFrame df {*chain};

    // Filter side events
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

    // Define distances in mm
    double base {256.};
    std::vector<double> dists;
    if(!isSide)
    {
        for(double d = 50; d < 110; d += 2)
            dists.push_back(base + d);
    }
    else if(mode == "l0")
    {
        for(double d = 60; d < 80; d += 2)
            dists.push_back(base + d);
    }
    else if(mode == "r0")
    {
        for(double d = 50; d < 100; d += 2)
            dists.push_back(0 - d);
    }

    int xbins {280};
    std::pair<double, double> xlims {-20, 400};
    int zbins {280};
    std::pair<double, double> zlims {-20, 400};
    // Save
    auto f {std::make_unique<TFile>(TString::Format("../Outputs/histos_%s.root", mode.Data()), "recreate")};
    std::cout << "Saving in : " << f->GetName() << '\n';
    f->WriteObject(&dists, "dists");
    TStopwatch timer {};
    timer.Start();
    for(const auto& dist : dists)
    {
        std::cout << "Distance : " << dist << '\n';
        // Redefine
        auto node {gated.Define("NewSP",
                                [&, dist](ActRoot::MergerData& d)
                                {
                                    // Reconstruct line from BP and SP
                                    auto p {d.fBP};
                                    auto dir {(d.fSP - d.fBP)};
                                    ActRoot::Line line {p, dir, 0};
                                    return isSide ? line.MoveToY(dist) : line.MoveToX(dist);
                                },
                                {"MergerData"})};
        // Fill histograms!
        ROOT::TThreadedObject<TH2D> hSP {ROOT::TNumSlots {node.GetNSlots()},
                                         "hSP",
                                         TString::Format("Side %.2f mm;Y [mm];Z [mm]", dist),
                                         xbins,
                                         xlims.first,
                                         xlims.second,
                                         zbins,
                                         zlims.first,
                                         zlims.second};
        std::map<int, ROOT::TThreadedObject<TH1D>> pxs, pzs;
        std::vector<int> idxs {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        for(const auto& idx : idxs)
        {
            pxs.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                        std::forward_as_tuple(ROOT::TNumSlots {node.GetNSlots()}, TString::Format("px%d", idx),
                                              TString::Format("%.2f mm X proj %d;X or Y [mm]", dist, idx), xbins,
                                              xlims.first, xlims.second));
            pzs.emplace(std::piecewise_construct, std::forward_as_tuple(idx),
                        std::forward_as_tuple(ROOT::TNumSlots {node.GetNSlots()}, TString::Format("pz%d", idx),
                                              TString::Format("%.2f mm Z proj %d;Z [mm]", dist, idx), zbins,
                                              zlims.first, zlims.second));
        }
        node.ForeachSlot(
            [&](unsigned int slot, ActRoot::MergerData& data, ROOT::Math::XYZPointF& sp)
            {
                hSP.GetAtSlot(slot)->Fill(isSide ? sp.X() : sp.Y(), sp.Z());
                auto idx {data.fSilNs.front()};
                if(pxs.count(idx))
                {
                    pxs[idx].GetAtSlot(slot)->Fill(isSide ? sp.X() : sp.Y());
                    pzs[idx].GetAtSlot(slot)->Fill(sp.Z());
                }
            },
            {"MergerData", "NewSP"});

        // Write data
        f->cd();
        auto path {TString::Format("d_%.1f_mm", dist)};
        auto* dir {f->mkdir(path)};
        dir->cd();
        hSP.Merge()->Write();
        for(auto m : {&pxs, &pzs})
        {
            for(auto& [_, h] : *m)
            {
                auto merged {h.Merge()};
                if(merged)
                    merged->Write();
            }
        }
    }


    timer.Stop();
    timer.Print();
}
