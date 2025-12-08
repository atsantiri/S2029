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


void understanding_pid()
{
    // Read data
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};

    ROOT::TThreadedObject<TH2D> hEsilThetaLab {
        "hEsilThetaLab", "ESil vs #theta_{Light}; #theta_{Light} [#circ]; E_{Sil} [MeV]", 200, 0, 180, 400, 0., 20};
    auto df = d.Filter([](ActRoot::MergerData& mer) { return mer.fLight.GetNLayers() == 1; }, {"MergerData"});

    df.Foreach([&](ActRoot::MergerData& m) { hEsilThetaLab->Fill(m.fThetaLight, m.fLight.fEs.front()); },
               {"MergerData"});

    auto* c1 {new TCanvas {"c1", "c1"}};
    hEsilThetaLab.Merge()->DrawClone("colz");
    auto l {new TLegend(0.7, 0.5, 0.9, 0.9)};

    for(int i = 0; i < 10; i++)
    {
        double ene = 65 - double(i) * 5;
        ActPhysics::Kinematics kin {TString::Format("17F(p,p)@%f", ene).Data()};
        auto* gp {kin.GetKinematicLine3()};
        gp->SetLineColor(i + 1);
        gp->Draw("same");
        l->AddEntry(gp, TString::Format("E_{beam} = %.0f MeV",ene));
    }
    l->Draw();
}