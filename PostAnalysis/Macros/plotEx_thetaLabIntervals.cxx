#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include "../HistConfig.h"
#include "../Settings.h"

// I want to plot the Ex for different theta lab values

void plotEx_thetaLabIntervals()
{
    // Read Pipe 3 output
    auto fIn {TString::Format("../Outputs/tree_ex_17F_p_p_%.2f_sil.root", S2029::EBeamIni)};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", fIn};

    // Create histograms
    auto hExTotal {df.Histo1D(HistConfig::Ex, "RecEx")};
    std::vector<ROOT::TThreadedObject<TH1D>*> hEx, hExR, hExL;
    auto hExModel {HistConfig::Ex};
    double step {10}; // deg
    double tmin {10};
    double tmax {70};
    int idx {0};
    for(double t = tmin; t < tmax; t += step)
    {
        hEx.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hEx%d", idx),
            TString::Format("#theta_{Lab} [%.1f, %.1f); E_{x} [MeV];Counts / %.f keV", t, t + step,
                            (hExModel.fXUp - hExModel.fXLow) / hExModel.fNbinsX * 1e3),
            hExModel.fNbinsX, hExModel.fXLow, hExModel.fXUp));
        hExL.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hEx%d", idx),
            TString::Format("#theta_{Lab} [%.1f, %.1f) L; E_{x} [MeV];Counts / %.f keV", t, t + step,
                            (hExModel.fXUp - hExModel.fXLow) / hExModel.fNbinsX * 1e3),
            hExModel.fNbinsX, hExModel.fXLow, hExModel.fXUp));
        hExR.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hEx%d", idx),
            TString::Format("#theta_{Lab} [%.1f, %.1f) R; E_{x} [MeV];Counts / %.f keV", t, t + step,
                            (hExModel.fXUp - hExModel.fXLow) / hExModel.fNbinsX * 1e3),
            hExModel.fNbinsX, hExModel.fXLow, hExModel.fXUp));
        idx++;
    }

    // Initialize slot 0 to not crash
    for(auto& h : hEx)
        h->GetAtSlot(0);

    // Fill histograms
    df.ForeachSlot(
        [&](unsigned int slot, const ActRoot::MergerData& m, double Ex)
        {
            for(int i = 0; i < hEx.size(); i++)
            {
                double t = tmin + i * step;
                if(m.fThetaLight >= t && m.fThetaLight < t + step)
                {
                    hEx[i]->GetAtSlot(slot)->Fill(Ex);
                    auto layer {m.fLight.fLayers.front()};
                    auto n {m.fLight.fNs.front()};
                    if(layer == "f0")
                    {
                        if(n == 0 || n == 1 || n == 4 || n == 5)
                            hExL[i]->GetAtSlot(slot)->Fill(Ex);

                        if(n == 2 || n == 3 || n == 6 || n == 7)
                            hExR[i]->GetAtSlot(slot)->Fill(Ex);
                    }
                    else if(layer == "l0")
                        hExL[i]->GetAtSlot(slot)->Fill(Ex);
                    else if(layer == "r0")
                        hExR[i]->GetAtSlot(slot)->Fill(Ex);
                }
            }
        },
        {"MergerData", "RecEx"});

    // Plot in canvas
    auto* c0 {new TCanvas("c0", "Ex for dthetaLab")};
    c0->DivideSquare(hEx.size());
    int c {1};

    auto* c1 {new TCanvas("c1", "Total Ex")};
    c1->cd();
    gStyle->SetPalette(91);
    hExTotal->DrawClone();
    auto* hs {new THStack};
    hs->SetTitle(TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV",
                                 (hExModel.fXUp - hExModel.fXLow) / hExModel.fNbinsX * 1e3));

    // for(auto& h : hEx)
    // {
    //     c0->cd(c);
    //     auto* clone {(TH1D*)h->Merge()->Clone()};
    //     clone->Draw();
    //     clone->SetLineWidth(2);
    //     hs->Add(clone);
    //     c++;
    // }

    for(int i = 0; i < idx; i++)
    {
        c0->cd(i + 1);
        auto* clone {(TH1D*)hEx[i]->Merge()->Clone()};
        clone->Draw();
        clone->SetLineWidth(2);
        hs->Add(clone);
        auto* cR {(TH1D*)hExR[i]->Merge()->Clone()};
        cR->SetLineColor(kMagenta);
        cR->Draw("same");
        cR->SetLineWidth(2);
        auto* cL {(TH1D*)hExL[i]->Merge()->Clone()};
        cL->SetLineColor(kOrange);
        cL->Draw("same");
        cL->SetLineWidth(2);
    }
    auto* leg1 = gPad->BuildLegend(0.6, 0.5, 0.9, 0.9);
    c1->cd();
    hs->Draw("plc pmc");
    auto* leg = gPad->BuildLegend(0.6, 0.5, 0.9, 0.9);
}