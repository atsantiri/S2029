#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"

#include "../PostAnalysis/HistConfig.h"
#include "../PostAnalysis/Settings.h"

// I want to test the drift velocity dependence in the reconstruction of Ex. By taking one thetalab look at phi~0 and
// the phi >0 and phi <0. The reconstructed Ex should in theory be the same.

void testEx_differentPhi()
{
    // Read Pipe 3 output
    auto fIn {TString::Format("../PostAnalysis/Outputs/tree_ex_17F_p_p_%.2f_sil.root",S2029::EBeamIni)};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {"Final_Tree", fIn};

    // angles
    double phimin {100.};
    double phimax {150.};
    double phi {50.};
    std::cout << "Phi: " << phi << std::endl;
    double thetamin {15.};
    double thetamax {25.};

    // Filter data
    auto df =
        d.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"; },
                 {"MergerData"}); // F0 data only
    auto dfCentre = df.Filter(
        [&](const ActRoot::MergerData& m)
        {
            if(m.fThetaLight >= thetamin && m.fThetaLight < thetamax)
            {
                if(TMath::Abs(m.fPhiLight) <= phi)
                    return true;
            }
            return false;
        },
        {"MergerData"});
    auto dfTop = df.Filter(
        [&](const ActRoot::MergerData& m)
        {
            if(m.fThetaLight >= thetamin && m.fThetaLight < thetamax)
            {
                if(m.fPhiLight > phimin && m.fPhiLight < phimax)
                    return true;
            }
            return false;
        },
        {"MergerData"});
    auto dfBot = df.Filter(
        [&](const ActRoot::MergerData& m)
        {
            if(m.fThetaLight >= thetamin && m.fThetaLight < thetamax)
            {
                if(m.fPhiLight < -phimin && m.fPhiLight > -phimax)
                    return true;
            }
            return false;
        },
        {"MergerData"});

    // Create histograms
    auto hExTotal {df.Filter([&](const ActRoot::MergerData& m)
                             { return (m.fThetaLight >= thetamin && m.fThetaLight < thetamax); }, {"MergerData"})
                       .Histo1D(HistConfig::Ex, "RecEx")};
    hExTotal->SetTitle(TString::Format("Ex for #theta_{Lab} [%.1f, %.1f)", thetamin, thetamax));
    auto hExCentre {dfCentre.Histo1D(HistConfig::Ex, "RecEx")};
    auto hExTop {dfTop.Histo1D(HistConfig::Ex, "RecEx")};
    auto hExBot {dfBot.Histo1D(HistConfig::Ex, "RecEx")};

    // auto hTheta {dfCentre.Histo1D({"hTheta", ";#theta_{Lab};Counts", 180, 0, 180}, "fThetaLight")};
    auto hPhi {df.Histo1D({"hPhi", ";#phi_{Lab};Counts", 200, -200, 200}, "fPhiLight")};
    auto hPhiRPx {df.Histo2D({"hPhiRPx", "#phi vs RPx for F0;RPx (mm);#phi_{Lab}", 280, 0, 280, 180, -180, 180},
                             "fRP.fCoordinates.fX", "fPhiLight")};
    auto hThetaRPx {df.Histo2D({"hThetaRPx", "#theta vs RPx;RPx (mm);#theta_{Lab}", 280, 0, 280, 180, 0, 180},
                             "fRP.fCoordinates.fX", "fThetaLight")};
    auto hPhiTheta {df.Histo2D({"hPhiTheta", "#phi vs #theta for F0 ;#theta_{Lab};#phi_{Lab}", 180, 0, 180, 180, -180, 180},
                             "fThetaLight", "fPhiLight")};

    auto* c0 {new TCanvas("c0", "Phi vs RPx")};
    hPhiRPx->DrawClone("colz");
    // hThetaRPx->DrawClone("colz");
    std::vector<double> angles = {phi, -phi, phimin, phimax, -phimin, -phimax};
    std::vector<int> colors = {616, 616, 800 - 3, 800 - 3, 101, 101};
    for(int i = 0; i < angles.size(); i++)
    {
        TLine* l = new TLine(0, angles[i], 280, angles[i]);
        l->SetLineStyle(2);
        l->SetLineColor(colors[i]);
        l->SetLineWidth(2);
        l->Draw("same");
    }

    auto* c1 {new TCanvas("c1", "Ex")};
    auto htot = (TH1*)hExTotal->Clone();
    htot->Draw();
    auto hc = (TH1*)hExCentre->Clone();
    hc->SetLineColor(kMagenta);
    hc->Draw("same");
    auto ht = (TH1*)hExTop->Clone();
    ht->SetLineColor(kOrange - 3);
    ht->Draw("same");
    auto hb = (TH1*)hExBot->Clone();
    hb->SetLineColor(kBrown);
    hb->Draw("same");
    auto leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(htot, "All", "l");
    leg->AddEntry(hc, "Centre", "l");
    leg->AddEntry(ht, "Top", "l");
    leg->AddEntry(hb, "Bot", "l");
    leg->Draw();
}