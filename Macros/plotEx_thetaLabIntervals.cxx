#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include "../PostAnalysis/HistConfig.h"

void plotEx_thetaLabIntervals()
{
    // Read Pipe 3 output
    auto fIn {TString("../PostAnalysis/Outputs/tree_ex_17F_p_p_3.90_sil.root")};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", fIn};

    // Create histograms
    auto hExTotal {df.Histo1D(HistConfig::Ex, "RecEx")};
    std::vector<ROOT::TThreadedObject<TH1D>*> hEx;
    auto hExModel {HistConfig::Ex};
    double step {10}; // deg
    double tmin {5};
    double tmax {85};
    int idx {0};
    for(double t = tmin; t < tmax; t += step)
    {
        hEx.push_back(new ROOT::TThreadedObject<TH1D>(
            TString::Format("hEx%d", idx),
            TString::Format("#theta_{Lab} [%.1f, %.1f); E_{x} [MeV];Counts / %.f keV", t, t + step,
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
                    hEx[i]->GetAtSlot(slot)->Fill(Ex);
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

    for(auto& h : hEx)
    {
        c0->cd(c);
        auto* clone {(TH1D*)h->Merge()->Clone()};
        clone->Draw();
        clone->SetLineWidth(2);
        hs->Add(clone);
        c++;
    }
    c1->cd();
    hs->Draw("plc pmc");
    auto* leg = gPad->BuildLegend(0.6, 0.5, 0.9, 0.9);
}