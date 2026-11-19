#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include "../HistConfig.h"
#include "../Settings.h"

//  make Ex(θlab) for all individual silicons

void debugF0()
{
    // Read Pipe 3 output
    auto fIn {TString::Format("../Outputs/tree_ex_17F_p_p_%.2f_sil.root", S2029::EBeamIni)};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"Final_Tree", fIn};

    // Create histograms
    int nsils {12};
    auto hExTotal {df.Histo1D(HistConfig::Ex, "RecEx")};
    std::map<int, std::map<int, std::unique_ptr<ROOT::TThreadedObject<TH1D>>>> hs;
    auto hExModel {HistConfig::Ex};
    double step {10}; // deg
    double tmin {5};
    double tmax {40};
    int idx {0};
    for(int s = 0; s < nsils; s++)
    {
        idx = 0;
        for(double t = tmin; t < tmax; t += step)
        {
            hs[s][idx] = std::make_unique<ROOT::TThreadedObject<TH1D>>(
                TString::Format("hEx_F0_%d_%d", s, idx),
                TString::Format("F0_%d, #theta_{Lab} [%.1f, %.1f); E_{x} [MeV];Counts / %.f keV", s, t, t + step,
                                (hExModel.fXUp - hExModel.fXLow) / hExModel.fNbinsX * 1e3),
                hExModel.fNbinsX, hExModel.fXLow, hExModel.fXUp);
            idx++;
        }
    }

    // Initialize slot 0 to not crash
    for(auto& [s, hEx] : hs)
    {
        for(auto& [idx, h] : hEx)
            h->GetAtSlot(0);
    }

    // Fill histograms
    df.ForeachSlot(
        [&](unsigned int slot, const ActRoot::MergerData& m, double Ex)
        {
            auto layer {m.fLight.fLayers.front()};
            auto n {m.fLight.fNs.front()};

            // only look at F0
            if(layer != "f0")
                return;

            for(int i = 0; i < hs[n].size(); i++)
            {
                double t = tmin + i * step;
                if(m.fThetaLight >= t && m.fThetaLight < t + step)
                {
                    hs[n][i]->GetAtSlot(slot)->Fill(Ex);
                    break;
                }
            }
        },
        {"MergerData", "RecEx"});

    // Plot in canvas
    std::vector<int> colors {TColor::GetColor("#4E79A7"), TColor::GetColor("#F28E2B"), TColor::GetColor("#E15759"),
                             TColor::GetColor("#76B7B2"), TColor::GetColor("#59A14F"), TColor::GetColor("#EDC948"),
                             TColor::GetColor("#B07AA1"), TColor::GetColor("#FF9DA7"), TColor::GetColor("#9C755F"),
                             TColor::GetColor("#BAB0AC"), TColor::GetColor("#499894"), TColor::GetColor("#D37295")};
    // for(int idx = 0; idx < hs.begin()->second.size(); idx++)
    // {
    //     c0->cd(c);
    //     int s {0};
    //     double maxY {0};
    //     std::vector<TH1D*> clones;
    //     for(auto& [det, hEx] : hs)
    //     {
    //         auto* clone {(TH1D*)hEx[idx]->Merge()->Clone(TString::Format("hEx_F0_%d_%d", det, idx))};

    //         clone->SetLineWidth(2);
    //         clone->Rebin(2);
    //         clone->SetLineColor(colors[det]);

    //         maxY = std::max(maxY, clone->GetMaximum());
    //         clones.push_back(clone);
    //     }
    //     for(auto* clone : clones)
    //     {
    //         clone->SetMaximum(maxY * 1.1);
    //         clone->Draw(s++ == 0 ? "" : "SAME");
    //     }

    //     c++;
    // }
    // auto* leg = gPad->BuildLegend(0.6, 0.5, 0.9, 0.9);
    std::map<int, std::map<int, TH1D*>> hMerged;

    for(auto& [det, hEx] : hs)
    {
        for(auto& [idx, h] : hEx)
        {
            hMerged[det][idx] = (TH1D*)h->Merge()->Clone(TString::Format("hEx_det%d_angle%d", det, idx));
        }
    }

    for(int idx = 0; idx < hs.begin()->second.size(); idx++)
    {
        auto* cDet {
            new TCanvas(TString::Format("cDet_%d", idx),
                        TString::Format("Ex, theta = [%.0f, %.0f) deg", tmin + idx * step, tmin + (idx + 1) * step))};

        cDet->DivideSquare(nsils);

        int pad {1};

        for(auto& [det, hEx] : hs)
        {
            cDet->cd(pad++);
            hMerged[det][idx]->SetLineWidth(2);
            hMerged[det][idx]->SetLineColor(colors[det]);
            hMerged[det][idx]->Draw();
        }
    }

    auto* c0 {new TCanvas("c0", "Ex for dthetaLab", 1200,1000)};
    c0->DivideSquare(hs.begin()->second.size());
    int c {1};

    for(int idx = 0; idx < hs.begin()->second.size(); idx++)
    {
        c0->cd(c);
        auto* stack {new THStack(TString::Format("stack_%d", idx),
                                 TString::Format("#theta_{Lab} [%.0f, %.0f); E_{x} [MeV]; Counts", tmin + idx * step,
                                                 tmin + (idx + 1) * step))};
        for(auto& [det, hEx] : hs)
        {
            hMerged[det][idx]->SetLineWidth(2);
            hMerged[det][idx]->SetLineColor(colors[det]);
            stack->Add(hMerged[det][idx]);
        }
        stack->Draw("hist");
        c++;
    }
    auto* leg = gPad->BuildLegend(0.6, 0.5, 0.9, 0.9);
}