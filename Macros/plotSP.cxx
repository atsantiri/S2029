#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActSilSpecs.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include <map>

void makeGrid(std::string layer)
{
    ActPhysics::SilSpecs specs;
    specs.ReadFile("../configs/silspecs.conf");
    auto sm {specs.GetLayer(layer).GetSilMatrix()->Clone()};
    auto& cuts = sm->GetGraphs();
    for(const auto& [id, cut] : cuts)
    {
        if(cut)
            cut->Draw("same");
    }
}

void plotSP()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    ROOT::RDataFrame df {*chain};

    // Fill histograms
    int nsils {12};
    std::map<std::string, std::map<int, ROOT::TThreadedObject<TH2D>>> hs;
    // Histogram model
    auto* h2d {new TH2D {"h2d", ";X or Y [mm];Z [mm]", 200, -10, 300, 200, -10, 300}};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        for(int s = 0; s < nsils; s++)
        {
            hs[layer].emplace(s, *h2d);
            hs[layer][s]->SetTitle(TString::Format("Layer %s", layer));
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
            auto sp {m.fLight.fSP};

            if(hs.count(layer))
            {
                if(hs[layer].count(n))
                {
                    if(layer == "f0")
                        hs[layer][n].Get()->Fill(sp.Y(), sp.Z());
                    else
                        hs[layer][n].Get()->Fill(sp.X(), sp.Z());
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
        // auto c = new TCanvas {lname, lname, 800, 600};
        TCanvas* c = nullptr;

        if(layer == "l0" || layer == "r0")
        {
            c = new TCanvas(lname, lname, 600, 800);
            c->Divide(3, 4);
        }
        else if(layer == "f0")
        {
            c = new TCanvas(lname, lname, 800, 600);
            c->Divide(4, 3);
        }
        // c0->cd(p);
        int idx {};
        std::cout << hsils.size() << " histograms for layer " << lname << std::endl;
        for(auto& [s, h] : hsils)
        {
            c->cd(idx + 1);
            auto color {idx + 1};
            if(color == 10) // 10 is white, as well as 0
                color = 46;
            auto opts {(idx == 0) ? "BOX" : "BOX same"};
            // Merge histos from threads
            h.Merge();
            // Set color
            h->SetMarkerColor(color);
            h->SetLineColor(color);
            // Set size
            h->SetMarkerStyle(6);
            // Set title
            h->SetTitle(TString::Format("%s_%d", lname.Data(), s));
            // Draw
            h->DrawClone("BOX");

            // make grid for comparison
            makeGrid(layer);
            idx++;
        }
    }
}
