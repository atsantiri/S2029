#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"

#include <map>
#include <string>

void simple_Pipe1_PID(const std::string& beam, const std::string& target, const std::string& light)
{
    std::string dataconf {};
    dataconf = "./../configs/data.conf";

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};

    // Only looking at L1 events
    auto df = d.Filter([](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                       { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); });

    ROOT::TThreadedObject<TH2D> hl1 {"hl1", "L1 PID;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1Gated {"hl1", "L1 PID > 100#circ;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0,
                                          3e5};
    ROOT::TThreadedObject<TH2D> hl1theta {
        "hl1theta", "L1 #theta;#theta_{L1} [#circ];Q_{total} [au]", 240, 0, 180, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1thetaCorr {
        "hl1thetaCorr", "L1 #thetas;#theta_{Light} [#circ];#theta_{Heavy} [#circ]", 240, 0, 180, 200, 0, 100};

    df.Foreach(
        [&](ActRoot::MergerData& m)
        {
            hl1->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
            hl1theta->Fill(m.fThetaLight, m.fLight.fQtotal);
            hl1thetaCorr->Fill(m.fThetaLight, m.fThetaHeavy);
            if(m.fThetaLight > 100)
                hl1Gated->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
        },
        {"MergerData"});
    auto* c0 {new TCanvas {"c0", "Pipe 1 PID L1"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hl1->DrawClone("colz");
    c0->cd(2);
    hl1theta->DrawClone("colz");
    c0->cd(3);
    hl1thetaCorr->DrawClone("colz");
    c0->cd(4);
    hl1Gated->DrawClone("colz");

}