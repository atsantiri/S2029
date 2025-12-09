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
    // auto filename {TString::Format("../Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(),
    // light.c_str())};


    std::ofstream streamer {"Outputs/eVertex_side_walls2.dat"};


    // std::string cutname {"pid_l1_what_is_this"};
    // auto filename {TString::Format("../Outputs/%s.root",cutname.c_str())};
    auto filename {TString::Format("../Outputs/tree_Ex_side_walls2.root")};
    std::cout << filename << std::endl;
    ROOT::RDataFrame df("Final_Tree", filename);
    // auto cols = df.GetColumnNames();
    // for(auto&& c : cols)
    //     std::cout << c << std::endl;

    // auto lambda {[](float x)
    //              {
    //                  double min {110}; // modify based on resonance location
    //                  double max {135};
    //                  return (min <= x) && (x <= max);
    //              }};


    // std::ofstream streamer {
    //     // TString::Format("./Outputs/events_%s_%s_%s.txt", beam.c_str(), target.c_str(), light.c_str())};
    //     TString::Format("./Outputs/%s.txt",cutname.c_str())};

    df.Foreach(
        [&](ActRoot::MergerData& mer, float x)
        {
            // if(lambda(x))
            mer.Stream(streamer);
        },
        {"MergerData", "fRP.fCoordinates.fX"});

    streamer.close();
}