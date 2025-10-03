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

void L1_Pipe1_PID()
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
    ROOT::RDataFrame df {*chain};

    // Only looking at L1 events
    auto lambdaL1 {[](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                   { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }};

    ROOT::TThreadedObject<TH2D> hl1 {"hl1", "L1 PID Light;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1thetaQ {"hl1thetaQ", "L1 #theta;#theta_{L1} [#circ];Q_{total} [au]", 240, 0, 180, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1thetaTL {"hl1thetaTL", "L1 #theta;#theta_{L1} [#circ];Raw TL [au]", 240, 0, 180, 200, 0, 120};
    ROOT::TThreadedObject<TH2D> hl1thetaCorr {"hl1thetaCorr", "L1 #thetas;#theta_{Light} [#circ];#theta_{Heavy} [#circ]", 240, 0, 180, 200, 0, 100};

    df.Foreach(
        [&](ActRoot::MergerData& m, ActRoot::ModularData& mod)
        {
            if(lambdaL1(m, mod))
            {
                hl1->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
                hl1thetaQ->Fill(m.fThetaLight, m.fLight.fQtotal);
                hl1thetaTL->Fill(m.fThetaLight, m.fLight.fRawTL);
                hl1thetaCorr->Fill(m.fThetaLight, m.fThetaHeavy);
            }
        },
        {"MergerData", "ModularData"});
    auto* c0 {new TCanvas {"c0", "Pipe 1 PID L1"}};
    c0->DivideSquare(4);
    c0->cd(1);
    hl1.Merge()->DrawClone("colz");
    c0->cd(2);
    hl1thetaQ.Merge()->DrawClone("colz");
    c0->cd(3);
    hl1thetaTL.Merge()->DrawClone("colz");
    c0->cd(4);
    hl1thetaCorr.Merge()->DrawClone("colz");


    ActRoot::CutsManager<std::string> cuts;
    std::string cutname {"pid_other_diagonal"};
    cuts.ReadCut("l1", TString::Format("./Cuts/%s.root", cutname.c_str()).Data());

    auto gated {df.Filter(
        [&](ActRoot::MergerData& m, ActRoot::ModularData& mod)
        {
            if(cuts.GetCut("l1"))
                return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal);
                // return cuts.IsInside("l1", m.fThetaLight, m.fThetaHeavy);
            else
                return false;
        },
        {"MergerData", "ModularData"})};
    auto fout {TString::Format("./Outputs/%s.root", cutname.c_str())};
    std::cout << "Saving PID_Tree in file : " << fout << '\n';
    gated.Snapshot("PID_Tree", fout.Data());

    c0->cd(1);
    cuts.DrawCut("l1");
}