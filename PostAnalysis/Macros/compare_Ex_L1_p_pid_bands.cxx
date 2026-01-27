#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"

#include "../HistConfig.h"


void compare_Ex_L1_p_pid_bands()
{
    TString ftop, fbot, fsil;
    ftop = TString::Format("../Outputs/tree_ex_17F_p_p_top_diag.root");
    fbot = TString::Format("../Outputs/tree_ex_17F_p_p_bot_diag.root");
    fsil = TString::Format("../Outputs/tree_ex_17F_p_p_3.84.root");


    ROOT::EnableImplicitMT();
    ROOT::RDataFrame dtop {"Final_Tree", ftop};
    ROOT::RDataFrame dbot {"Final_Tree", fbot};
    ROOT::RDataFrame dsil {"Final_Tree", fsil};

    auto hExTop {dtop.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "RecEx")};

    auto hExBot {dbot.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "RecEx")};

    auto hExSil {dsil.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "RecEx")};

    auto* c {new TCanvas("c", "", 1000, 800)};
    c->SetLogy();
    hExBot->SetLineColor(2);
    hExSil->SetLineColor(4);
    hExBot->DrawClone();
    hExTop->DrawClone("same");
    hExSil->DrawClone("same");
}