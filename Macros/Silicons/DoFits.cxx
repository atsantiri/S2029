#include "ActSilMatrix.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>

#include "./GetContourFuncs.cxx"

void DoFits(TString mode)
{
    mode.ToLower();
    std::cout << "Running for : " << mode << '\n';

    // Assign things depending on mode
    std::vector<int> idxs {0, 1, 2, 3, 4, 5, 6, 7, 9, 11};
    bool isSide {mode.Contains("f") ? false : true};
    std::cout << "isSide ? " << std::boolalpha << isSide << '\n';
    if(mode == "l0")
        idxs = {1, 2, 4, 5, 7, 8, 10, 11};
    if(mode == "r0")
        idxs = {0, 1, 4, 6, 7, 9, 10};

    // Read data
    auto fin {std::make_unique<TFile>(TString::Format("../Outputs/%s_histograms.root", mode.Data()).Data())};
    // Set histograms to process
    std::map<int, TH1D*> pxys, pzs;
    for(auto& idx : idxs)
    {
        auto xykey {TString::Format("pxy%d", idx)};
        auto zkey {TString::Format("pz%d", idx)};
        pxys[idx] = fin->Get<TH1D>(xykey);
        pxys[idx]->SetDirectory(nullptr);
        pzs[idx] = fin->Get<TH1D>(zkey);
        pzs[idx]->SetDirectory(nullptr);
    }

    ProjMap nxys, nzs;
    for(const auto& idx : idxs)
    {
        double w {5};
        double s {1};
        double scalew {15};
        // XY
        auto* fx {FindBestFit(pxys[idx], w, s)};
        nxys[idx] = ScaleWithFunc(pxys[idx], fx, scalew);
        // Z
        auto* fz {FindBestFit(pzs[idx], w, s)};
        nzs[idx] = ScaleWithFunc(pzs[idx], fz, scalew);
    }


    // Fit to contour
    double thresh {0.85};
    double width {10};
    if(isSide)
    {
        thresh = 0.6;
        width = 20; // increase for lateral sils
    }
    // XY
    auto xypoints {FitToCountour(nxys, thresh, width)};
    // Z
    auto zpoints {FitToCountour(nzs, thresh, width)};

    // Build graphs with resolutions
    auto* gsy {new TGraphErrors};
    gsy->SetTitle("#sigma_{XY};Sil idx;#sigma_{XY} [mm]");
    auto* gsz {new TGraphErrors};
    gsz->SetTitle("#sigma_{Z};Sil idx;#sigma_{Z} [mm]");
    // Y
    for(auto& [idx, h] : nxys)
    {
        // Get functions
        auto* left {h->GetFunction("fleft")};
        auto* right {h->GetFunction("fright")};
        if(left && right)
        {
            auto sum {left->GetParError(2) + right->GetParameter(2)};
            gsy->AddPoint(idx, sum / 2);
        }
    }
    // Z
    for(auto& [idx, h] : nzs)
    {
        // Get functions
        auto* left {h->GetFunction("fleft")};
        auto* right {h->GetFunction("fright")};
        if(left && right)
        {
            auto sum {left->GetParError(2) + right->GetParameter(2)};
            gsz->AddPoint(idx, sum / 2);
        }
    }
    // Fit them
    for(auto* g : {gsy, gsz})
    {
        g->Fit("pol0", "0Q");
        g->GetFunction("pol0")->ResetBit(TF1::kNotDraw);
    }


    // Build class to store them
    auto* smatrix {new ActPhysics::SilMatrix {mode.Data()}};
    for(auto& [i, ypoint] : xypoints)
    {
        const auto& zpoint {zpoints[i]};
        smatrix->AddSil(i, ypoint, zpoint);
    }

    // plot
    auto* cpx {new TCanvas("cpx", "XY projection canvas")};
    PlotAll(cpx, pxys);

    auto* cnx {new TCanvas("cnx", "Normalized XY canvas")};
    PlotAll(cnx, nxys);

    auto* cpz {new TCanvas("cpz", "Z projection canvas")};
    PlotAll(cpz, pzs);

    auto* cnz {new TCanvas("cnz", "Normalized Z canvas")};
    PlotAll(cnz, nzs);

    gStyle->SetOptFit();
    auto* cs {new TCanvas {"cs", "Sigma canvas"}};
    cs->DivideSquare(2);
    for(auto* g : {gsy, gsz})
    {
        g->SetMarkerStyle(24);
        g->SetLineWidth(2);
    }
    cs->cd(1);
    gsy->Draw("apl");
    cs->cd(2);
    gsz->Draw("apl");

    auto* cm {new TCanvas("cm", "Sil matrix")};
    auto name {TString::Format("../Outputs/%s_matrix.root", mode.Data())};
    std::cout << "Saving matrix at : " << name << '\n';
    smatrix->Write(name.Data());
    smatrix->Draw(false);
}
