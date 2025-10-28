#include "ActDataManager.h"
#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"

#include <map>

void plotAllSP()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    ROOT::RDataFrame df {*chain};

    // Fill histograms
    int nsils {12};
    std::map<std::string, std::map<int, ROOT::TThreadedObject<TH2D>>> hs;
    // Histogram model
    auto* h2d {new TH2D {"h2d", ";X or Y [pad];Z [btb]", 300, 0, 300, 500, 0, 500}};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        for(int s = 0; s < nsils; s++)
        {
            hs[layer].emplace(s, *h2d);
        }
    }
    int rstat {0}, lstat {0};
    df.Foreach(
        [&](ActRoot::MergerData& m)
        {
            if(m.fLight.GetNLayers() != 1)
                return;
            // No need to check for SP, it has already been done in Filter
            auto layer {m.fLight.fLayers.front()}; // ensured to have at least size >= 1
            auto n {m.fLight.fNs.front()};
            auto sp {m.fLight.fSP};
            if(hs.count(layer))
            {
                if(hs[layer].count(n))
                {
                    if(layer == "f0")
                        hs[layer][n].Get()->Fill(sp.Y(), sp.Z());
                    else
                    {
                        hs[layer][n].Get()->Fill(sp.X(), sp.Z());
                        (layer == "l0") ? lstat++
                                        : rstat++; // keep track of stats for L and R to check they're symmetric
                    }
                }
            }
        },
        {"MergerData"});

    //// Draw
    int canvasIdx {0};
    for(auto& [layer, hsils] : hs)
    {
        auto lname {TString(layer)};
        lname.ToUpper();
        auto c = new TCanvas {lname, lname, 800, 600};

        int idx {0};
        std::cout << hsils.size() << " histograms for layer " << lname << std::endl;
        for(auto& [s, h] : hsils)
        {
            auto color {idx + 1};
            if(color == 10) // 10 is white, as well as 0
                color = 46;
            auto opts {(idx == 0) ? "BOX" : "BOX same"};
            // Merge histos from threads
            h.Merge();
            h->SetStats(false);
            h->SetTitle(lname);
            // Set color
            h->SetMarkerColor(color);
            h->SetLineColor(color);
            // Set size
            h->SetMarkerStyle(6);
            // Draw
            h->DrawClone(opts);
            idx++;
        }
    }
    std::cout << "Number of counts for R0: " << rstat << " and L0: " << lstat << std::endl;
}
