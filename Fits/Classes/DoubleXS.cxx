// Function from https://github.com/loopset/S2008/tree/main/Fits 
#ifndef DoubleXS_cxx
#define DoubleXS_cxx
#include "DoubleXS.h"

#include "ActKinematics.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"

#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>

DoubleXS::DoubleXS(TH2* hData, TH2* hEff, ActPhysics::SRIM* srim, double nb, double rho, ActPhysics::Kinematics* kin,
                   TString isCM)
    : fHist((TH2*)hData->Clone("fHist")),
      fEff(hEff),
      fsrim(srim),
      fNbeams(nb),
      fDensity(rho),
      fKin(kin),
      fIsCM(isCM)
{
    fOriginal = (TH2*)hData->Clone("hOriginal");
    fHist->SetTitle("Cross section");
    fHist->GetZaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");
    fHist->GetZaxis()->SetMaxDigits(2);
    fHist->Sumw2(); // store sum of squares of weights
    ApplyEff();
    ApplySolidAngle();
    ApplyThickness();
    ApplyNormalisation();
}

void DoubleXS::ApplyEff()
{
    // Build TProfile
    fEffProf = new TProfile2D {
        "hEffProf",
        TString::Format("Mean efficiency;#theta_{%s} [#circ];E_{%s} [MeV];#epsilon", fIsCM.Data(), fIsCM.Data()),
        fHist->GetNbinsX(),
        fHist->GetXaxis()->GetXmin(),
        fHist->GetXaxis()->GetXmax(),
        fHist->GetNbinsY(),
        fHist->GetYaxis()->GetXmin(),
        fHist->GetYaxis()->GetXmax()};
    // Fill it
    for(int x = 1; x <= fEff->GetNbinsX(); x++)
    {
        for(int y = 1; y <= fEff->GetNbinsY(); y++)
        {
            auto content {fEff->GetBinContent(x, y)};
            if(content <= 0)
                continue;
            auto angle {fEff->GetXaxis()->GetBinCenter(x)};
            auto e {fEff->GetYaxis()->GetBinCenter(y)};
            fEffProf->Fill(angle, e, content);
        }
    }
    // Divide
    fHist->Divide(fEffProf);
}

void DoubleXS::ApplySolidAngle()
{

    auto binOmega {fHist->GetXaxis()->GetBinUpEdge(1)}; // in deg
    fFuncOmega = new TF1 {
        "fFuncOmega",
        TString::Format("TMath::TwoPi() * TMath::Sin(x * TMath::DegToRad()) * (%f) * TMath::DegToRad()", binOmega),
        fHist->GetXaxis()->GetXmin(), fHist->GetXaxis()->GetXmax()};
    fFuncOmega->SetTitle(
        TString::Format("Solid angle func;#theta_{%s} [#circ];2#pi sin(#theta_{%s}) #Delta #theta_{%s} [sr]",
                        fIsCM.Data(), fIsCM.Data(), fIsCM.Data()));
    fHist->Divide(fFuncOmega);
}

double DoubleXS::TransformCMtoLab(double ecm)
{
    auto m1 {fKin->GetParticle(1).GetMass()};
    auto m2 {fKin->GetParticle(2).GetMass()};
    return (m1 + m2) / m2 * ecm;
}

void DoubleXS::ApplyThickness()
{
    fThick = new TH2D {"Thick",
                       TString::Format("SRIM thickness per E_{%s} bin;#theta_{%s} [#circ];E_{%s} [MeV];Thick [cm]",
                                       fIsCM.Data(), fIsCM.Data(), fIsCM.Data()),
                       fHist->GetNbinsX(),
                       fHist->GetXaxis()->GetXmin(),
                       fHist->GetXaxis()->GetXmax(),
                       fHist->GetNbinsY(),
                       fHist->GetYaxis()->GetXmin(),
                       fHist->GetYaxis()->GetXmax()};
    if(!fsrim)
        return;
    for(int y = 1; y <= fThick->GetNbinsY(); y++)
    {
        // This is in direct kinematics
        auto elow {fThick->GetYaxis()->GetBinLowEdge(y)};
        auto eup {fThick->GetYaxis()->GetBinUpEdge(y)};
        // Get equivalent 20Mg energy as beam
        if(!fKin)
            throw std::runtime_error("Cannot transform p energy to RIB");
        elow = fKin->ComputeEquivalentOtherT1(elow);
        eup = fKin->ComputeEquivalentOtherT1(eup);
        // Convert to lab if in CM
        if(fIsCM == "CM")
        {
            elow = TransformCMtoLab(elow);
            eup = TransformCMtoLab(eup);
        }
        auto key {"beam"};
        if(!fsrim || !fsrim->CheckKeyIsStored(key))
            throw std::runtime_error("Cannot access SRIM pointer");

        auto rlow {fsrim->EvalRange(key, elow)};
        auto rup {fsrim->EvalRange(key, eup)};
        // And DeltaX
        auto deltax {TMath::Abs(rup - rlow)};

        // Fill that value for all X bins
        for(int x = 1; x <= fThick->GetNbinsX(); x++)
        {
            fThick->SetBinContent(x, y, deltax / 10); // mm to cm (later on, density is in atoms/cm3)
            fThick->SetBinError(x, y, 0);             // this has no error (in principle...)
        }
    }
    fHist->Divide(fThick);
}

void DoubleXS::ApplyNormalisation()
{
    // Divide by density
    fHist->Scale(1. / fDensity);
    // Divide by Nbeams
    fHist->Scale(1. / fNbeams);
    // Convert to [mb/sr]
    fHist->Scale(1e27);
}

void DoubleXS::Draw()
{
    auto* c0 {new TCanvas {"cxsCalc", "Double XS calculations", 1000,800}};
    c0->DivideSquare(6);
    c0->cd(1);
    fOriginal->Draw("colz");
    c0->cd(2);
    fEff->Draw("colz");
    c0->cd(3);
    fEffProf->Draw("colz1");
    c0->cd(4);
    fFuncOmega->Draw();
    c0->cd(5);
    fThick->Draw("colz");
    c0->cd(6);
    fHist->Draw("colz");
}

void DoubleXS::Project(int cThresh)
{
    // Project onto thetaCM axis
    int idx {};
    for(int y = 1; y <= fHist->GetNbinsY(); y++)
    {
        auto proj {fHist->ProjectionX(TString::Format("pTheta%d", idx), y, y)};
        if(proj->GetEntries() < cThresh)
        {
            delete proj;
            continue;
        }
        else
        {
            auto low {fHist->GetYaxis()->GetBinLowEdge(y)};
            auto up {fHist->GetYaxis()->GetBinUpEdge(y)};
            proj->SetTitle(TString::Format("E_{%s} #in [%.2f,%.2f) MeV;#theta_{%s} [#circ];d#sigma/d#Omega [mb/sr]",
                                           fIsCM.Data(), low, up, fIsCM.Data()));
            fProjsThetaCM.push_back(proj);
            fIvsThetaCM.push_back({low, up});
            idx++;
        }
    }
    // Project onto ECM axis
    idx = 0;
    for(int x = 1; x <= fHist->GetNbinsX(); x++)
    {
        auto proj {fHist->ProjectionY(TString::Format("pE%d", idx), x, x)};
        if(proj->GetEntries() < cThresh)
        {
            delete proj;
            continue;
        }
        else
        {
            auto low {fHist->GetXaxis()->GetBinLowEdge(x)};
            auto up {fHist->GetXaxis()->GetBinUpEdge(x)};
            proj->SetTitle(TString::Format("#theta_{%s} #in [%.2f,%.2f)#circ;E_{%s} [MeV];d#sigma/d#Omega [mb/sr]",
                                           fIsCM.Data(), low, up, fIsCM.Data()));
            fProjsECM.push_back(proj);
            fIvsECM.push_back({low, up});
            idx++;
        }
    }
}

void DoubleXS::DrawProjectionsThetaCM(const std::function<void(TH1*)>& apply)
{
    auto* c {new TCanvas {"cxsTheta", "Theta projection canvas"}};
    c->DivideSquare(fProjsThetaCM.size());
    for(int i = 0; i < fProjsThetaCM.size(); i++)
    {
        c->cd(i + 1);
        auto& h {fProjsThetaCM[i]};
        if(apply)
            apply(h);
        h->Draw("histe");
    }
}

void DoubleXS::DrawProjectionsECM(const std::function<void(TH1*)>& apply)
{
    auto* c {new TCanvas {"cxsE", "E projection canvas",1500,1000}};
    c->DivideSquare(fProjsECM.size());
    for(int i = 0; i < fProjsECM.size(); i++)
    {
        c->cd(i + 1);
        auto& h {fProjsECM[i]};
        if(apply)
            apply(h);
        h->Draw("histe");
    }
}

TH1D* DoubleXS::GetProjectionECM(double thetamin, double thetamax)
{
    auto bmin {fHist->GetXaxis()->FindBin(thetamin)};
    auto min {fHist->GetXaxis()->GetBinCenter(bmin)};
    auto bmax {fHist->GetXaxis()->FindBin(thetamax) - 1};
    auto max {fHist->GetXaxis()->GetBinCenter(bmax)};
    auto* p {fHist->ProfileY("pInRange", bmin, bmax)};
    p->SetTitle(TString::Format("#theta_{%s} #in [%.2f,%.2f)", fIsCM.Data(), min, max));
    return p;
}

void DoubleXS::WriteInAzureFormat(int idx, const TString& file, TH1D* pout, const std::pair<double, double>& ivs)
{
    // Ensure projection exists
    TH1D* p {};
    if(pout)
        p = pout;
    else
    {
        try
        {
            p = fProjsECM.at(idx);
        }
        catch(std::exception& e)
        {
            throw std::runtime_error("DoubleXS::WriteInAzureFormat: cannot find idx: " + std::to_string(idx) +
                                     " in fProjecECM");
        }
    }
    double centre {};
    if(pout)
        centre = (ivs.first + ivs.second) / 2;
    else
        centre = (fIvsECM[idx].first + fIvsECM[idx].second) / 2;
    std::ofstream streamer {file};
    for(int b = 1; b <= p->GetNbinsX(); b++)
    {
        auto x {p->GetBinCenter(b)};
        auto y {p->GetBinContent(b)};
        if(y <= 0)
            continue;
        auto uy {p->GetBinError(b)};
        // Workaround for azure...
        if(y - uy <= 0)
            uy = y - 0.01;
        streamer << x << "    " << centre << "    " << y << "    " << uy << '\n';
    }
    streamer.close();
}

#endif // DoubleXS_cxx
