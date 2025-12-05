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

    auto* c1 {new TCanvas {"c1", "Theta Lab vs CM elastic"}};
    auto* gThetaLabvsCM {kp.GetThetaLabvsThetaCMLine()};
    gThetaLabvsCM->Draw();

    auto* c2 {new TCanvas {"c2", "Particle 3 vs Particle 4"}};
    c2->DivideSquare(3);
    // c2->cd(1);
    // auto* gTheoTheta3vsTheta4 {kp.GetTheta3vs4Line()};
    // gTheoTheta3vsTheta4->Draw();
    // c2->cd(2);
    // gp->Draw();
    // c2->cd(3);
    // auto* gp4 {kp.GetKinematicLine4()};
    // gp4->SetLineColor(2);
    // gp4->Draw();

    for(int i = 0; i < 10; i++)
    {
        double ene = 65.5 - double(i) * 2.;
        std::cout<<ene<<std::endl;
        ActPhysics::Kinematics kin {TString::Format("17F(p,p)@%f", ene).Data()};
        c2->cd(1);
        auto* gTheoTheta3vsTheta4 {kin.GetTheta3vs4Line()};
        gTheoTheta3vsTheta4->SetLineColor(i + 1);
        if(i == 0)
            gTheoTheta3vsTheta4->Draw();
        else
            gTheoTheta3vsTheta4->Draw("same");

        c2->cd(2);
        auto* gp {kin.GetKinematicLine3()};
        gp->SetLineColor(i + 1);
        if(i == 0)
            gp->Draw();
        else
            gp->Draw("same");
        c2->cd(3);
        auto* gp4 {kin.GetKinematicLine4()};
        gp4->SetLineColor(i + 1);
        if(i == 0)
            gp4->Draw();
        else
            gp4->Draw("same");
    }
}
