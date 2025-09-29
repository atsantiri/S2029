#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include <fstream>


void extractEvtNumber()
{
    std::string beam {"17F"};
    std::string target {"p"};
    std::string light {"4He"};
    auto filename {TString::Format("../Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};

    ROOT::RDataFrame df("PID_Tree", filename);
    // auto cols = df.GetColumnNames();
    // for(auto&& c : cols)
    //     std::cout << c << std::endl;


    auto gateRPx {df.Filter(
        [](float x)
        {
            double min {110}; // modify based on resonance location
            double max {135};
            return (min <= x) && (x <= max);
        },
        {"fRP.fCoordinates.fX"})};

    std::ofstream streamer {
        TString::Format("./Outputs/events_%s_%s_%s.txt", beam.c_str(), target.c_str(), light.c_str())};
    gateRPx.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();
}