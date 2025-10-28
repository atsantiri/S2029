#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"


int silSpectra()
{
    std::string dataconf {};
    dataconf = "./../configs/data.conf";

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};

    // Filter Si hit
    auto df {d.Filter([](ActRoot::MergerData& mer) { return mer.fLight.GetNLayers() == 1; }, {"MergerData"})};

    int nsils {12};
    std::map<std::string, std::map<int, ROOT::TThreadedObject<TH1D>>> hs;

    auto hSil {new TH1D {"hSil", ";E_{Sil} [MeV];Counts", 200, 0, 20}};

    for(const auto& layer : {"f0", "l0", "r0"})
    {
        for(int s = 0; s < nsils; s++)
        {
            hs[layer].emplace(s, *hSil);
            hs[layer][s]->SetTitle(TString::Format("%s_%d", layer, s));
        }
    }

    df.Foreach(
        [&](ActRoot::MergerData& m)
        {
            if(m.fLight.GetNLayers() != 1)
                return;
            // No need to check for SP, it has already been done in Filter
            auto layer {m.fLight.fLayers.front()}; // ensured to have at least size >= 1
            auto n {m.fLight.fNs.front()};

            if(hs.count(layer) && hs[layer].count(n))
                hs[layer][n].Get()->Fill(m.fLight.fEs.front());
        },
        {"MergerData"});

int canvasIdx {0};
    for(auto& [layer, hsils] : hs)
    {
        auto c = new TCanvas {layer.c_str(), layer.c_str(), 800, 600};
        if(layer == "l0" || layer == "r0")
            c->Divide(3, 4);
        if(layer == "f0")
            c->Divide(4, 3);
        int idx {};
        for(auto& [s, h] : hsils)
        {
            c->cd(idx + 1);
            // Merge histos from threads
            h.Merge();
            // Set title
            h->SetTitle(TString::Format("%s_%d", layer.c_str(), s));
            // Draw
            h->DrawClone();
            idx++;
        }
    }    

    return 0;
}