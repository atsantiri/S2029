#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <fstream>
#include <stdexcept>

void count_L1ok_events()
{
    std::string dataconf {};
    dataconf = "./../configs/data.conf";

    // ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman {dataconf, ActRoot::ModeType::EReadSilMod};
    auto chain {datman.GetJoinedData()};
    auto chain2 {datman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());

    ROOT::RDataFrame df {*chain};

    // Get GATCONF values
    auto def {df.Define("GATCONF", [](ActRoot::ModularData& mod) { return static_cast<int>(mod.fLeaves["GATCONF"]); },
                        {"ModularData"})
                  .Define("INCONF", [](ActRoot::ModularData& mod) { return static_cast<int>(mod.fLeaves["INCONF"]); },
                          {"ModularData"})};

    // count triggers
    std::atomic<unsigned long int> l1ok_gat {};
    std::atomic<unsigned long int> l1ok_inc {};

    std::ofstream streamer {"./Outputs/events_L1ok.txt"};

    def.Foreach(
        [&](const int& gatconf, const int& inconf, ActRoot::MergerData& mer)
        {
            if(gatconf == 8)
            {
                l1ok_gat++;
                mer.Stream(streamer);
            }
            if(inconf == (64 + 8) && (gatconf == 64))
                l1ok_inc++;
        },
        {"GATCONF", "INCONF", "MergerData"});


    // Print report
    std::cout << "===== GATCONF report =====" << '\n';
    std::cout << "-> L1Ok from GATCONF only " << l1ok_gat << '\n';
    std::cout << "-> L1Ok from INCONF  " << l1ok_inc << '\n';
    std::cout << "==========================" << '\n';
    std::cout << "L1Ok event list written in ./Outputs/events_L1ok.txt" << std::endl;
}
