#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"

void beamCum()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    ROOT::RDataFrame df {*chain};

    auto gated {df.Filter([&](ActRoot::TPCData& tpc, ActRoot::ModularData& mod)
                          { return (mod.Get("GATCONF") == 64) && (tpc.fClusters.size() == 1); },
                          {"TPCData", "ModularData"})};

    auto* h2d {new TH2D {"hPad", "Pad plane;X;Y", 128, 0, 128, 128, 0, 128}};
    gated.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            auto& cluster {tpc.fClusters.front()};
            for(const auto& v : cluster.GetVoxels())
            {
                auto& pos {v.GetPosition()};
                h2d->Fill(pos.X(), pos.Y());
            }
        },
        {"TPCData"});

    auto* c0 {new TCanvas {"c0", "beam cum"}};
    h2d->Draw("colz");
}
