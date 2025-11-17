#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TROOT.h"
#include "TString.h"

#include <fstream>

void findFunEvts()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());

    ROOT::RDataFrame df {*chain};

    auto df_filtered {
        df.Filter(
            [](ActRoot::TPCData& tpc, ActRoot::ModularData& m)
            {
                return (tpc.fRPs.size() == 1 && tpc.fClusters.size() > 4);
            }, 
            {"TPCData", "ModularData"})
        // .Filter(
        //     [](ActRoot::TPCData& tpc)
        //     {
        //         auto& rp {tpc.fRPs[0]};
        //         double min {110}; // modify based on resonance location
        //         double max {140};
        //         return (min <= rp.X()) && (rp.X() <= max);
        // },
        // {"TPCData"})
    };

    std::ofstream streamer {"../Outputs/funEvts_many_clusters.txt"};
    df_filtered.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();

}