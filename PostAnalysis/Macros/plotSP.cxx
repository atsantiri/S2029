#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"

#include <map>

void plotSP()
{
    ROOT::RDataFrame df {"Final_Tree", "../Outputs/tree_ex_17F_p_p.root"};

    // Fill histograms
    int nsils {12};
    std::map<std::string, std::map<int, ROOT::TThreadedObject<TH2D>>> hs;
    // Histogram model
    auto* h2d {new TH2D {"h2d", "SP;X or Y [pad];Z [btb]", 300, 0, 300, 500, 0, 500}};
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

    // Draw
    auto* c0 {new TCanvas {"c0", "SP canvas"}};
    c0->DivideSquare(4);
    int p {1};
    int canvasIdx {0};
    for(auto& [layer, hsils] : hs)
    {
        // Crear un nuevo canvas para cada histograma
        auto cname = Form("c%d", canvasIdx++);
        auto c = new TCanvas {cname, Form("SP canvas %d", canvasIdx), 800, 600};
        if(layer == "l0" || layer == "r0")
            c->Divide(3, 4);
        if(layer == "f0")
            c->Divide(4, 3);
        // c0->cd(p);
        int idx {};
        std::cout << hsils.size() << " histograms for layer " << layer << std::endl;
        for(auto& [s, h] : hsils)
        {
            c->cd(idx + 1);
            auto color {idx + 1};
            if(color == 10) // 10 is white, as well as 0
                color = 46;
            auto opts {(idx == 0) ? "scat" : "scat same"};
            // Merge histos from threads
            h.Merge();
            // Set color
            h.GetAtSlot(0)->SetMarkerColor(color);
            // Set size
            h.GetAtSlot(0)->SetMarkerStyle(6);
            // Set title
            h->SetTitle(TString::Format("Layer %s, S%d", layer.c_str(), s));
            // Draw
            h.GetAtSlot(0)->DrawClone("scat");
            idx++;
        }
        p++;
    }
}
