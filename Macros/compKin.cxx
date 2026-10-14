#include "ActKinematics.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include "../PostAnalysis/Settings.h"
void compKin()
{
    ActPhysics::Particle pb {"17F"};
    auto Ebeam = S2029::EBeamIni * pb.GetAMU();
    ActPhysics::Kinematics kp {TString::Format("17F(p,p)@%.2f",Ebeam).Data()};
    ActPhysics::Kinematics ka {TString::Format("17F(p,a)@%.2f",Ebeam).Data()};

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
    c1->DivideSquare(2);
    c1->cd(1);
    auto* gThetaLabvsCM {kp.GetThetaLabvsThetaCMLine()};
    gThetaLabvsCM->GetYaxis()->SetRangeUser(0,90);
    c1->cd(1)->SetGridx();
    c1->cd(1)->SetGridy();
    gThetaLabvsCM->Draw("AL");


    auto* c2 {new TCanvas {"c2", "Particle 3 vs Particle 4"}};
    c2->DivideSquare(2);

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
        double ene = Ebeam - double(i) * 2.5;
        std::cout << ene << std::endl;
        ActPhysics::Kinematics kin {TString::Format("17F(p,p)@%f", ene).Data()};
        c2->cd(1);
        auto* gTheoTheta3vsTheta4 {kin.GetTheta3vs4Line()};
        gTheoTheta3vsTheta4->SetLineColor(i + 1);
        if(i == 0)
            gTheoTheta3vsTheta4->Draw("apl plc pmc");
        else
            gTheoTheta3vsTheta4->Draw("apl plc pmc same");

        c2->cd(2);
        auto* gp4 {kin.GetKinematicLine4()};
        gp4->SetLineColor(i + 1);
        if(i == 0)
            gp4->Draw();
        else
            gp4->Draw("same");

        c1->cd(2);
        auto* gp {kin.GetKinematicLine3()};
        gp->SetLineColor(i + 1);
        if(i == 0)
            gp->Draw();
        else
            gp->Draw("same");
    }
}
