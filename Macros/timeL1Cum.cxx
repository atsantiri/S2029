#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"
#include "ActTPCParameters.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"

void timeL1Cum()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    ActRoot::TPCParameters fTPC {"Actar"};
    int rebinZ {fTPC->GetREBINZ()};

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};
    auto df {d.Filter(
        [](ActRoot::MergerData& mer, ActRoot::ModularData& mod, ActRoot::TPCData& tpc)
        {
            return ((mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1)) ||
                   ((mod.Get("GATCONF") == 64) && (tpc.fClusters.size() == 1));
        },
        {"MergerData", "ModularData", "TPCData"})};

    ROOT::TThreadedObject<TH2D> h2d {"hPad", "Pad plane;X [pad];Y [pad]", 128, 0, 128, 128, 0, 128};
    ROOT::TThreadedObject<TH2D> hside {"hSide", "Side;X [pad];Z [tb]", 128, 0, 128, 512, 0, 512};
    df.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            if(tpc.fClusters.empty())
                return;
            auto& cluster {tpc.fClusters.front()};
            for(const auto& v : cluster.GetVoxels())
            {
                auto& pos {v.GetPosition()};
                h2d->Fill(pos.X(), pos.Y());
                hside->Fill(pos.X(), pos.Z());
            }
        },
        {"TPCData"});

    auto* c0 {new TCanvas {"c0", "L1 events"}};
    c0->DivideSquare(2);
    c0->cd(1);
    h2d.Merge()->DrawClone("colz");
    c0->cd(2);
    hside.Merge()->DrawClone("colz");
}
