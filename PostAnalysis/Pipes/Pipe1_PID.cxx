#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TString.h"

#include <map>
#include <string>

void Pipe1_PID(const std::string& beam, const std::string& target, const std::string& light)
{
    std::string dataconf {};
    dataconf = "./../configs/data.conf";

    // Read data
    ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    // RDataFrame
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // LIGHT particle
    // Define lambda functions
    // 1-> Stopped in first silicon layer
    auto lambdaOne {[](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }};
    // 2-> In two layers
    auto lambdaTwo {[](ActRoot::MergerData& m)
                    {
                        if(m.fLight.GetNLayers() == 2)
                            return (m.fLight.GetLayer(0) == "f0" && m.fLight.GetLayer(1) == "f1");
                        else
                            return false;
                    }};
    // 3-> L1
    auto lambdaIsL1 {[](ActRoot::MergerData& mer, ActRoot::ModularData& mod)
                     { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }};

    // Fill histograms
    std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgas, hstwo, hszero;
    // Histogram models
    auto hGasSil {new TH2D {"hGasSil", ";E_{Sil} [MeV];#Delta E_{gas} [arb. units]", 450, 0, 70, 600, 0, 3000}};
    auto hTwoSils {new TH2D {"hTwoSils", ";#DeltaE_{0} [MeV];#DeltaE_{1} [MeV]", 500, 0, 80, 400, 0, 30}};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        hsgas.emplace(layer, *hGasSil);
        hsgas[layer]->SetTitle(TString::Format("%s", layer));
    }
    hstwo.emplace("f0-f1", *hTwoSils);
    hstwo["f0-f1"]->SetTitle("f0-f1");
    ROOT::TThreadedObject<TH2D> hl1 {"hl1", "L1 PID;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH1D> hCHIO {"hCHIO", "CHIO;Channel [adc];Counts [au]", 1000, 0, 16384};
    ROOT::TThreadedObject<TH2D> hTACCFA {"hTACCFA", "CFA TAC; CHIO Channel [adc]; TAC [adc]", 1000,0,16384, 1000, 0, 16384};

    ROOT::TThreadedObject<TH2D> hl1Gated {"hl1", "L1 PID > 100#circ;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0,
                                          3e5};
    ROOT::TThreadedObject<TH2D> hl1theta {
        "hl1theta", "L1 #theta;#theta_{L1} [#circ];Q_{total} [au]", 240, 0, 180, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1thetaCorr {
        "hl1thetaCorr", "L1 #thetas;#theta_{Light} [#circ];#theta_{Heavy} [#circ]", 240, 0, 180, 200, 0, 100};

    // Fill them
    df.Foreach(
        [&](ActRoot::MergerData& m, ActRoot::ModularData& mod)
        {
            // L1
            if(lambdaIsL1(m, mod))
            {
                hl1->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
                hl1theta->Fill(m.fThetaLight, m.fLight.fQtotal);
                hl1thetaCorr->Fill(m.fThetaLight, m.fThetaHeavy);
                if(m.fThetaLight > 100)
                    hl1Gated->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
                return;
            }
            // Light
            if(lambdaOne(m)) // Gas-E0 PID
            {
                auto layer {m.fLight.GetLayer(0)};
                if(hsgas.count(layer))
                    hsgas[layer]->Fill(m.fLight.fEs.front(), m.fLight.fQave);
            }
            else if(lambdaTwo(m)) // E0-E1 PID
            {
                hstwo["f0-f1"]->Fill(m.fLight.fEs[0], m.fLight.fEs[1]);
            }

            hCHIO->Fill(mod.Get("CHIO"));
            hTACCFA->Fill(mod.fLeaves["CHIO"],mod.fLeaves["TAC_CFA_HF"]);
            
            
        },
        {"MergerData", "ModularData"});

    // If cuts are present, apply them
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("l0", TString::Format("./Cuts/pid_%s_l0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("r0", TString::Format("./Cuts/pid_%s_r0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("f0", TString::Format("./Cuts/pid_%s_f0_%s.root", light.c_str(), beam.c_str()).Data());
    cuts.ReadCut("l1", TString::Format("./Cuts/pid_%s_l1_%s.root", light.c_str(), beam.c_str()).Data());
    
    
    // Two sils PID
    // cuts.ReadCut("f0-f1", TString::Format("./Cuts/pid_%s_f0_f1_%s.root", light.c_str(), beam.c_str()).Data());
    // Get list of cuts
    auto listOfCuts {cuts.GetListOfKeys()};
    if(listOfCuts.size())
    {
        // Apply PID and save in file
        auto gated {df.Filter(
            [&](ActRoot::MergerData& m, ActRoot::ModularData& mod)
            {
                // L1
                if(lambdaIsL1(m, mod))
                {
                    if(cuts.GetCut("l1"))
                        return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal);
                    else
                        return false;
                }
                // One silicon
                else if(lambdaOne(m))
                {
                    auto layer {m.fLight.GetLayer(0)};
                    if(cuts.GetCut(layer))
                    {
                        // LIGHT particle
                        auto l {cuts.IsInside(layer, m.fLight.fEs[0], m.fLight.fQave)};
                        
                        return l ;
                    }
                    else
                        return false;
                }
                else if(cuts.GetCut("f0-f1") && lambdaTwo(m)) // PID in fo-f1
                    return cuts.IsInside("f0-f1", m.fLight.fEs[0], m.fLight.fEs[1]);
                else
                    return false;
            },
            {"MergerData", "ModularData"})};
        auto name {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
        std::cout << "Saving PID_Tree in file : " << name << '\n';
        gated.Snapshot("PID_Tree", name.Data());
    }

    // Draw
    auto* c0 {new TCanvas {"c0", "Pipe 1 PID canvas 0"}};
    c0->DivideSquare(6);
    int p {1};
    c0->cd(1);
    for(auto& [layer, h] : hsgas)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz"); //merge histograms written by all threads
        cuts.DrawCut(layer);
        p++;
    }
    for(auto& [layer, h] : hstwo)
    {
        c0->cd(p);
        h.Merge()->DrawClone("colz");
        p++;
    }
    c0->cd(5);
    // hl1.Merge()->DrawClone("colz");
    // cuts.DrawCut("l1");
    // c0->cd(6);
    // hl1theta.Merge()->DrawClone("colz");

    // auto* c2 {new TCanvas {"c2", "Pipe1 PID canvas 2"}};
    // c2->DivideSquare(4);
    // c2->cd(1);
    // hl1.GetAtSlot(0)->DrawClone("colz");
    // c2->cd(2);
    // hl1theta.GetAtSlot(0)->DrawClone("colz");
    // c2->cd(3);
    // hl1thetaCorr.Merge()->DrawClone("colz");
    // // auto* gtheo {ActPhysics::Kinematics(TString::Format("%s(d,p)@82.5", beam.c_str()).Data()).GetTheta3vs4Line()};
    // gtheo->SetLineColor(46);
    // gtheo->Draw("l");
    // c2->cd(4);
    // hl1Gated.Merge()->DrawClone("colz");



    auto* c3 {new TCanvas {"c3", "Pipe1 IC and CFA TAC"}};
    c3->DivideSquare(2);
    c3->cd(1);
    hCHIO.Merge()->DrawClone();
    c3->cd(2);
    hTACCFA.Merge()->DrawClone("COLZ");
}
