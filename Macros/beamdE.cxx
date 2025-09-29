//////////////////////////////////////////////////////////////////////////////////////////
// Script to investigate dE/dx from 17O and 17F beams from the two runs where the pads  //
// were polarized near the entrance of ACTAR                                            //
//////////////////////////////////////////////////////////////////////////////////////////


#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"

#include <fstream>

TH1D* make_dE(int run, TCanvas* c, TLegend* l, int pad)
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    dataman.SetRuns(run, run);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    ROOT::RDataFrame df {*chain};


    auto gated {df.Filter(
                      [](ActRoot::TPCData& tpc, ActRoot::ModularData& mod, ActRoot::MergerData& mer)
                      {
                          return (mod.Get("GATCONF") == 64) && (tpc.fClusters.size() == 1) && (mer.fEntry > 500) &&
                                 (mer.fEntry < 1354);
                      },
                      {"TPCData", "ModularData", "MergerData"})
                    .Define("GATCONF", [](ActRoot::ModularData& mod)
                            { return static_cast<int>(mod.fLeaves["GATCONF"]); }, {"ModularData"})};

    std::atomic<unsigned long int> cfa {};
    gated.Foreach(
        [&](const int& gatconf)
        {
            // if(gatconf == 64)
            cfa++;
        },
        {"GATCONF"});


    auto h2d = new TH2D {
        TString::Format("hPad%d", run), TString::Format("Pad plane %d;X [pad];Y [pad]", run), 128, 0, 128, 128, 0, 128};
    auto hdE = new TH1D {TString::Format("hdE%d", run), "; X [pad]; Q per pad", 30, 0, 30};
    gated.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            auto& cluster {tpc.fClusters.front()};
            for(const auto& v : cluster.GetVoxels())
            {
                auto posX {v.GetPosition().X()};
                auto posY {v.GetPosition().Y()};
                auto q {v.GetCharge()};
                h2d->Fill(posX, posY);

                if(posX <= 30)
                    hdE->Fill(posX, q);
            }
        },
        {"TPCData"});


    c->cd(pad);
    h2d->DrawClone("colz");
    hdE->SetLineColor(pad);
    std::cout << "Normalizing run " << run << " to number of incoming particles: " << cfa << std::endl;
    hdE->Scale(1. / cfa);
    l->AddEntry(hdE, Form("Run %d", run));
    c->cd(3);
    hdE->DrawClone("same");

    return hdE;
}

TH1D* do_simu(const std::string p, TCanvas* c, TLegend* l, int color)
{
    auto* srim {new ActPhysics::SRIM()};
    ActPhysics::Particle beam {p};

    auto load_table = [&](const char* label, const TString& path)
    {
        if(std::ifstream(path.Data()))
        {
            srim->ReadTable(label, path.Data());
        }
        else
        {
            std::cerr << path << " not found!\n";
            return false;
        }
        return true;
    };

    std::vector<std::pair<const char*, TString>> tables = {
        {"beamInMylar", TString::Format("../Simulation/SRIM/%s_mylar.txt", p.c_str())},
        {"beamIniC4H10", TString::Format("../Simulation/SRIM/%s_iC4H10_6mbar.txt", p.c_str())},
        {"beamInACTARgas", TString::Format("../Simulation/SRIM/%s_H2-iC4H10_95-5_700mbar.txt", p.c_str())}};

    for(const auto& [label, path] : tables)
    {
        if(!load_table(label, path))
            return nullptr;
    }

    int iter {1000};
    auto hSim = new TH1D {TString::Format("hSim%s", p.c_str()), "", 30, 0, 60};

    for(int i = 0; i < iter; i++)
    {
        auto Tini {5.5 * beam.GetAMU()};
        auto TCFAm {srim->Slow("beamInMylar", Tini, 4.8e-3)};
        auto TCFA {srim->Slow("beamIniC4H10", TCFAm, 9.6)};
        auto TAfterFoil {srim->Slow("beamInMylar", TCFA, 12e-3)};
        auto TActiveArea {srim->Slow("beamInACTARgas", TAfterFoil, 31.6 + 12.8 + 7.6)};
        double currentEne = TActiveArea;
        for(int x = 0; x < 60; x++)
        {
            double dEdx = srim->EvalStoppingPower("beamInACTARgas", currentEne);
            currentEne = srim->Slow("beamInACTARgas", currentEne, 1);
            hSim->Fill(x, dEdx);
        }
    }
    hSim->SetLineColor(color);
    hSim->Scale(1. / iter);
    hSim->Draw("same");
    l->AddEntry(hSim, p.c_str());

    return hSim;
}


void beamdE()
{

    bool compare_with_simu {false};
    // ROOT::EnableImplicitMT();
    auto* c {new TCanvas {"c", "Data", 1500, 800}};
    auto l = new TLegend(0.5, 0.1, 0.9, 0.3);

    c->Divide(3);
    auto htemp = new TH2D {"htemp", "; X [pad]; Q per pad", 100, 0, 30, 100, 500, 7000};
    htemp->SetStats(false);
    c->cd(3);
    htemp->Draw();

    auto hO = make_dE(43, c, l, 1); // 43 is 50V polarization before changing charged state
    auto hF = make_dE(45, c, l, 2); // 45 is with 9+ charged state at 70V polarization
    l->Draw();

    // How does the simulation with srim compare for the two?
    auto* c1 {new TCanvas {"c1", "Simu", 1000, 800}};
    c1->Divide(2);
    c1->cd(1);
    auto htemp2 = new TH2D {"htemp2", "; X [mm]; dE/dx", 100, 0, 60, 100, 0, .5};
    auto l1 = new TLegend(0.7, 0.8, 0.9, 0.9);
    htemp2->SetStats(false);
    htemp2->Draw();
    auto hsimO = do_simu("17O", c1, l1, 1);
    auto hsimF = do_simu("17F", c1, l1, 2);
    l1->Draw();


    // Find what the ratio b/w 17O and 17F should have been
    TH1D* hRatioSim = (TH1D*)hsimF->Clone("hRatioSim");
    hRatioSim->Divide(hsimO);
    c1->cd(2);
    hRatioSim->GetXaxis()->SetLimits(0, 30);
    hRatioSim->Draw();

    c->cd(3);
    TH1D* hOscaled = (TH1D*)hO->Clone("hOscaled");
    hOscaled->Multiply(hRatioSim);
    hOscaled->SetLineColor(3);
    hOscaled->Draw("same");
    l->AddEntry(hOscaled, "where F should have been");
    l->Draw();

    // Now what is the ratio of the corrected 17O tp 17F. Is it the 70V/50V of the different polarization in the two
    // runs?
    TH1D* hRatio2 = (TH1D*)hOscaled->Clone("hRatio2");
    hRatio2->Divide(hF);

    auto* c2 {new TCanvas {"c2", "", 800, 800}};
    c2->cd();
    hRatio2->SetLineColor(1);
    hRatio2->GetYaxis()->SetRangeUser(0, 2);
    hRatio2->Draw();
    // TH1D* hRatio2Corrected = (TH1D*)hRatio2->Clone("hRatio2Corrected");
    // hRatio2Corrected->Divide(hRatioSim);
    // hRatio2Corrected->SetLineColor(2);
    // hRatio2Corrected->Draw("same");

    auto l2 = new TLegend(0.7, 0.8, 0.9, 0.9);
    l2->AddEntry(hRatio2, "7/5?");
    TLine* line = new TLine(0, 7. / 5., 30, 7. / 5.);
    line->SetLineColor(kRed);
    line->Draw("same");
    l2->AddEntry(line, "7/5");
    l2->Draw();
}
