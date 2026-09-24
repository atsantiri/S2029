#include "ActSilMatrix.h"
#include "ActSilSpecs.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TPaveText.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <vector>

// #include "../../PostAnalysis/Utils.cxx"

void OverlapSM(const std::string& beam)
{
    std::cout << "Beam : " << beam << '\n';
    // Read histograms
    auto file {new TFile {TString::Format("./Outputs/histos_%s.root", beam.c_str())}};
    auto* hxz {file->Get<TH2D>("hTrajXZ")};
    auto* hyz {file->Get<TH2D>("hTrajYZ")};
    // Mean of beam in histogram
    auto meanFront {hyz->GetMean(2)};
    auto meanSide {hxz->GetMean(2)};

    // Real silicon specs
    ActPhysics::SilSpecs specs;
    specs.ReadFile("../../configs/silspecs.conf");

    // Map things
    std::vector<TH2D*> hs {hyz, hxz, hxz};
    // std::vector<ActPhysics::SilMatrix*> sms {S2008::GetFrontMatrix(), S2008::GetLeftMatrix(), S2008::GetRightMatrix()};
    std::vector<std::string> labels {"f0", "l0", "r0"};
    std::vector<std::string> layers {"f0", "l0", "r0"};
    std::vector<ActPhysics::SilMatrix*> phys;
    // Format phys sms
    for(int i = 0; i < labels.size(); i++)
    {
        auto sm {specs.GetLayer(layers[i]).GetSilMatrix()};
        phys.push_back(sm->Clone());
        phys.back()->SetName(labels[i]);
    }

    // Get means of desired silicons
    std::vector<double> zmeans;
    std::vector<double> diffs;
    std::vector<TPaveText*> texts;
    int idx {};
    for(auto* sm : phys)
    {
        TString label {labels[idx]};
        label.ToLower();
        std::vector<double> temp;
        std::set<int> sils;
        double ref {};
        if(label.Contains("f0"))
        {
            sils = {4, 7};
            ref = meanFront;
        }
        else if(label.Contains("l0"))
        {
            sils = {4};
            ref = meanSide;
        }
        else if(label.Contains("r0"))
        {
            sils = {4};
            ref = meanFront;
        }
        for(auto sil : sils)
        {
            double x {};
            double y {};
            sm->GetSil(sil)->Center(x, y);
            temp.push_back(y);
        }
        // Compute mean
        zmeans.push_back(std::accumulate(temp.begin(), temp.end(), 0.0) / temp.size());
        // And diff
        diffs.push_back(zmeans.back() - ref);
        auto* text {new TPaveText {0.4, 0.75, 0.6, 0.88, "NDC"}};
        text->AddText(TString::Format("#DeltaZ = %.2f mm", diffs.back()));
        text->SetBorderSize(0);
        texts.push_back(text);
        // Print:
        std::cout << "Mean for " << labels[idx] << " : " << zmeans.back() << " mm" << '\n';
        // Move center of physical sms to this value
        phys[idx]->MoveZTo(zmeans.back(), sils);
        idx++;
    }


    // Draw
    auto* c0 {new TCanvas {"c0", "SM and Emittance canvas"}};
    c0->DivideSquare(4);
    for(int i = 0; i < labels.size(); i++)
    {
        c0->cd(i + 1);
        hs[i]->DrawCopy()->SetTitle(labels[i].c_str());
        sms[i]->Draw();
        texts[i]->Draw();
    }

    auto* c1 {new TCanvas {"c1", "Physical silicons"}};
    c1->DivideSquare(4);
    for(int i = 0; i < phys.size(); i++)
    {
        c1->cd(i + 1);
        phys[i]->Draw(false);
        auto* cl {sms[i]->Clone()};
        cl->SetSyle(false, 0, 0, 3001);
        cl->Draw();
    }
}
