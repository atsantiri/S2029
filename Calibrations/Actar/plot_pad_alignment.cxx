#include "ActDataManager.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"

#include <vector>
#include "TCanvas.h"
#include "TProfile.h"


void plot_pad_alignment()
{
    auto* c0 {new TCanvas {"c0", "Pad alignment"}};
    c0->DivideSquare(5);

    int p {1};

    for (int run=12; run<=16; ++run) {
        auto df{ROOT::RDataFrame("GETTree", TString::Format("../../../RootFiles/Cluster/Clusters_Run_00%d.root",run))};

        // auto* hq {new TH1D {"hq", "pulser run;pad ID ; Qave", 16384, 0, 16383}};
        auto hq {new TProfile(TString::Format("hq%d", run), TString::Format("Pulser run %d;pad ID;Average Q(ch)",run), 16384, 0, 16383)}; //TProfile computes avrages instead of sum

        auto df2 {df.Define("padIDs", 
                    [](const ActRoot::TPCData& tpc) {
                        std::vector<int> pads;
                        for (const auto& v : tpc.fRaw)
                            pads.push_back((int)v.GetPosition().X() * 128 + (int)v.GetPosition().Y());
                        return pads;
                    }, 
                    {"TPCData"})
                .Define("charges", 
                    [](const ActRoot::TPCData& tpc) {
                        std::vector<double> q;
                        for (const auto& v : tpc.fRaw)
                            q.push_back(v.GetCharge());
                        return q;
                    }, 
                    {"TPCData"})};

        df2.Foreach(
            [&](const std::vector<int>& pads, const std::vector<double>& charges){
                for (size_t i=0; i<pads.size(); ++i)
                    hq->Fill(pads[i], charges[i]);  // TProfile stores average per bin
            },
            {"padIDs","charges"}
        );

        c0->cd(p++);
        // hq->SetLineColor(kOrange);
        hq->Draw();
        hq->GetYaxis()->SetRangeUser(0, 100);

        c0->cd(6);
        if (run==12){
            auto* hqAll = (TProfile*)hq->Clone("hqAll");
            hqAll->SetTitle("All runs");
            hqAll->Draw();
            }
        else hq->Draw("same");

    }

}
