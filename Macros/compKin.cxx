#include "ActKinematics.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"

void compKin()
{
    ActPhysics::Kinematics kp {"17F(p,p)@65.5"};
    ActPhysics::Kinematics ka {"17F(p,a)@65.5"};

    auto* mg {new TMultiGraph};
    mg->SetTitle("Kin comparison;#theta_{Lab} [#circ];E_{Lab} [MeV]");
    auto* gp {kp.GetKinematicLine3()};
    gp->SetTitle("17F(p,p)");
    auto* ga {ka.GetKinematicLine3()};
    ga->SetTitle("17F(p,a)");
    mg->Add(gp);
    mg->Add(ga);

    auto* c0 {new TCanvas {"c0", "Comp kin"}};
    mg->Draw("apl plc pmc");
    auto l = new TLegend(0.7, 0.7, 0.9, 0.9);
    l->AddEntry(gp, gp->GetTitle());
    l->AddEntry(ga, ga->GetTitle());
    l->Draw();


    auto* c1 {new TCanvas{"c1", "Theta Lab vs CM elastic"}};
    auto* gThetaLabvsCM {kp.GetThetaLabvsThetaCMLine()};
    gThetaLabvsCM->Draw();

}
