#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <fstream>

void p2p()
{
    // the three body vertices exist only in the filter data. So we start from there.
    ActRoot::DataManager dataman {"./configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    // The run numbers and event entry exist in the merger data so we want that too.
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    // We make the two chains friends so we can access the run and entry numbers from the filter data.
    chain->AddFriend(chain2.get());

    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // Gate on candidates for (p,2p)
    auto gated {df.Filter([](ActRoot::TPCData& tpc) { return (tpc.fClusters.size() == 4) && (tpc.fRPs.size() == 1); }, // we want the cluster 
                          {"TPC"})};
    // And also gate on X values
    auto gateRPx {gated.Filter(
        [](ActRoot::TPCData& tpc)
        {
            auto& rp {tpc.fRPs[0]};
            double min {45};// modify based on resonance location
            double max {90};
            return (min <= rp.X()) && (rp.X() <= max);
        },
        {"TPCData"})};

    // Plot histogram of RPx for event with 3 vertices. Is there some structure here we could use??
    auto hRPx {gated.Define("RPx", [](ActRoot::TPCData& tpc) { return tpc.fRPs[0].X(); }, {"TPCData"})
                   .Histo1D({"hRPx", "RP x histo;RP.X [mm or pads, check];Counts", 300, 0, 256}, "RPx")};

    // Write to file: run and entry numbers that passed the cuts for 3-body vertices within min and max RPx
    std::ofstream streamer {"./p2p_list.txt"};
    gateRPx.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();

    // Draw
    auto* c0 {new TCanvas {"c0", "(p,2p) canvas"}};
    hRPx->DrawClone();
}
