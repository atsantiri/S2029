#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"

#include <fstream>

void calcZOffset()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    ActPhysics::SilSpecs specs;
    specs.ReadFile("../../configs/silspecs.conf");

    ROOT::RDataFrame d {*chain};

    auto df {d.Filter([&](ActRoot::TPCData& tpc, ActRoot::ModularData& mod)
                      { return (mod.Get("GATCONF") == 64) && (tpc.fClusters.size() == 1); },
                      {"TPCData", "ModularData"})};

    std::vector<double> zpos;
    // auto* h2d {new TH2D {"hPad", "Pad plane;X [pad];Z [pad]", 128, 0, 128, 128, 0, 128}};
    df.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            auto& cluster {tpc.fClusters.front()};
            for(const auto& v : cluster.GetVoxels())
            {
                auto& pos {v.GetPosition()};
                zpos.push_back(pos.Z());
                // h2d->Fill(pos.X(), pos.Z());
            }
        },
        {"TPCData"});

    double mean_zpos {TMath::Mean(zpos.size(), zpos.data())};
    std::cout << "average zpos: " << mean_zpos << std::endl;

    // get fDriftFactor to convert pos.Z to mm
    ActRoot::InputParser parser {"../../configs/detector.conf"};
    auto block {parser.GetBlock("Merger")};
    auto fDriftFactor {block->GetDouble("DriftFactor")};
    auto mean_zpos_mm {mean_zpos*fDriftFactor};

    // auto* c0 {new TCanvas {"c0", "beam XZ"}};
    // h2d->Draw("colz");

    for(const auto& layer : {"f0", "l0", "r0"})
    {
        auto z {specs.GetLayer(layer).GetPoint().Z()};
        std::cout << layer << ": " << z << ", zOffser = "<< mean_zpos_mm - z<<std::endl;
    }
}