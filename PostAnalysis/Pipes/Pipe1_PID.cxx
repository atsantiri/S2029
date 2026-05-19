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
#include "TLegend.h"
#include "TString.h"

#include <map>
#include <string>

void Pipe1_PID(const std::string &beam, const std::string &target, const std::string &light)
{
    std::string dataconf{};
    dataconf = "./../configs/data.conf";

    // Read data
    // ActRoot::DataManager dataman {dataconf, ActRoot::ModeType::EMerge};
    // auto chain {dataman.GetChain()};
    // auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    // chain->AddFriend(chain2.get());

    // // RDataFrame
    // ROOT::EnableImplicitMT();
    // ROOT::RDataFrame df {*chain};
    auto filename{TString::Format("./Outputs/tree_preProcessed_%s.root", beam.c_str())};
    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df{"PreProcessed_Tree", filename};

    // LIGHT particle
    // Define lambda functions
    // 1-> Stopped in first silicon layer
    auto lambdaOne{[](ActRoot::MergerData &m)
                   { return m.fLight.GetNLayers() == 1; }};
    // 2-> In two layers
    auto lambdaTwo{[](ActRoot::MergerData &m)
                   {
                       if (m.fLight.GetNLayers() == 2)
                           return (m.fLight.GetLayer(0) == "f0" && m.fLight.GetLayer(1) == "f1");
                       else
                           return false;
                   }};
    // 3-> L1
    auto lambdaIsL1{[](ActRoot::MergerData &mer, ActRoot::ModularData &mod)
                    { return (mod.Get("GATCONF") == 8) && (mer.fLightIdx != -1); }};

    // Fill histograms
    std::map<std::string, ROOT::TThreadedObject<TH2D>> hsgas, hstwo, hszero;
    // Histogram models
    auto hGasSil{new TH2D{"hGasSil", ";E_{Sil} [MeV];#Delta E_{gas} [arb. units]", 400, 0, 40, 400, 0, 2000}};
    auto hTwoSils{new TH2D{"hTwoSils", ";#DeltaE_{0} [MeV];#DeltaE_{1} [MeV]", 500, 0, 80, 400, 0, 30}};
    for (const auto &layer : {"f0", "l0", "r0"})
    {
        hsgas.emplace(layer, *hGasSil);
        hsgas[layer]->SetTitle(TString::Format("%s", layer));
    }
    hstwo.emplace("f0-f1", *hTwoSils);
    hstwo["f0-f1"]->SetTitle("f0-f1");
    ROOT::TThreadedObject<TH2D> hEsilThetaLab{
        "hEsilThetaLab", "ESil vs #theta_{Light}; #theta_{Light} [#circ]; E_{Sil} [MeV]", 200, 0, 180, 400, 0., 20};

    ROOT::TThreadedObject<TH2D> hl1{"hl1", "L1 PID;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1Gated{"hl1", "L1 PID > 100#circ;Raw TL [au];Q_{total} [au]", 200, 0, 120, 2000, 0,
                                         3e5};
    ROOT::TThreadedObject<TH2D> hl1theta{
        "hl1theta", "L1 #theta;#theta_{L1} [#circ];Q_{total} [au]", 240, 0, 180, 2000, 0, 3e5};
    ROOT::TThreadedObject<TH2D> hl1thetaCorr{
        "hl1thetaCorr", "L1 #thetas;#theta_{Light} [#circ];#theta_{Heavy} [#circ]", 240, 0, 180, 200, 0, 100};

    // Fill them
    df.Foreach(
        [&](ActRoot::MergerData &m, ActRoot::ModularData &mod)
        {
            // L1
            if (lambdaIsL1(m, mod))
            {
                // if(m.fPhiLight < 2)
                // {
                hl1->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
                hl1theta->Fill(m.fThetaLight, m.fLight.fQtotal);
                hl1thetaCorr->Fill(m.fThetaLight, m.fThetaHeavy);
                if (m.fThetaLight > 100)
                    hl1Gated->Fill(m.fLight.fRawTL, m.fLight.fQtotal);
                // }
                return;
            }
            // Light
            if (lambdaOne(m)) // Gas-E0 PID
            {
                auto layer{m.fLight.GetLayer(0)};
                if (hsgas.count(layer))
                {
                    hsgas[layer]->Fill(m.fLight.fEs.front(), m.fLight.fQave);
                }
                hEsilThetaLab->Fill(m.fThetaLight, m.fLight.fEs.front());
            }
            else if (lambdaTwo(m)) // E0-E1 PID
            {
                hstwo["f0-f1"]->Fill(m.fLight.fEs[0], m.fLight.fEs[1]);
            }
        },
        {"MergerData", "ModularData"});

    // If cuts are present, apply them
    ActRoot::CutsManager<std::string> cuts;

    // Gas PID
    // cuts.ReadCut("l0", TString::Format("./Cuts/pid_%s_l0_%s.root", light.c_str(), beam.c_str()).Data());
    // cuts.ReadCut("r0", TString::Format("./Cuts/pid_%s_r0_%s.root", light.c_str(), beam.c_str()).Data());
    // cuts.ReadCut("f0", TString::Format("./Cuts/pid_%s_f0_%s.root", light.c_str(), beam.c_str()).Data());
    // cuts.ReadCut("l1", TString::Format("./Cuts/pid_%s_l1_%s.root", light.c_str(), beam.c_str()).Data());
    // std::string cutname {"pid_l1_veryLowQevents"};
    std::string cutname{"pid_4He_l1_17F"};
    // std::string cutname{"pid_l1_lowQ_alphas"};
    cuts.ReadCut("l1", TString::Format("./Cuts/%s.root", cutname.c_str()).Data());
    cuts.ReadCut("l1p", "./Cuts/pid_p_l1_17F.root");
    cuts.ReadCut("l1bad", "./Cuts/pid_lowQ_region.root");

    // Two sils PID
    // cuts.ReadCut("f0-f1", TString::Format("./Cuts/pid_%s_f0_f1_%s.root", light.c_str(), beam.c_str()).Data());
    // Get list of cuts
    auto listOfCuts{cuts.GetListOfKeys()};
    if (listOfCuts.size())
    {
        // Apply PID and save in file
        auto gated{df.Filter(
            [&](ActRoot::MergerData &m, ActRoot::ModularData &mod)
            {
                // L1
                if (lambdaIsL1(m, mod))
                {
                    if (cuts.GetCut("l1"))
                        return cuts.IsInside("l1", m.fLight.fRawTL, m.fLight.fQtotal);

                    else
                        return false;
                }
                // One silicon
                else if (lambdaOne(m))
                {
                    auto layer{m.fLight.GetLayer(0)};
                    if (cuts.GetCut(layer))
                    {
                        // LIGHT particle
                        auto l{cuts.IsInside(layer, m.fLight.fEs[0], m.fLight.fQave)};

                        return l;
                    }
                    else
                        return false;
                }
                else if (cuts.GetCut("f0-f1") && lambdaTwo(m)) // PID in fo-f1
                    return cuts.IsInside("f0-f1", m.fLight.fEs[0], m.fLight.fEs[1]);
                else
                    return false;
            },
            {"MergerData", "ModularData"})};
        // auto name {TString::Format("./Outputs/grass.root")};
        // auto name {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
        auto name{TString::Format("./Outputs/tree_%s.root", cutname.c_str())};
        std::cout << "Saving PID_Tree in file : " << name << '\n';
        gated.Snapshot("PID_Tree", name.Data());
        auto tot = gated.Count();
        std::cout << "Counts in cut " << cutname << " " << tot.GetValue() << std::endl;
    }

    std::map<std::string, ROOT::RDF::RResultPtr<ULong64_t>> cutCounts;
    for (const auto &cut : cuts.GetListOfKeys())
    {
        cutCounts[cut] = df.Filter(
                               [&, cut](ActRoot::MergerData &m, ActRoot::ModularData &mod)
                               {
                                   if (cut == "l1" || cut == "l1p" || cut == "l1bad")
                                   {
                                       return lambdaIsL1(m, mod) &&
                                              cuts.IsInside(cut,
                                                            m.fLight.fRawTL,
                                                            m.fLight.fQtotal);
                                   }
                                   return false;
                               },
                               {"MergerData", "ModularData"})
                             .Count();
    }
    std::cout << "\n==== COUNTS PER CUT ====\n";
    for (auto &[cut, counter] : cutCounts)
        std::cout << cut << " : " << counter.GetValue() << '\n';

    // Draw
    auto *c11{new TCanvas{"c11", "Pipe 1 PID Si"}};
    c11->DivideSquare(4);
    int p{1};
    c11->cd(1);
    for (auto &[layer, h] : hsgas)
    {
        c11->cd(p);
        h.Merge()->DrawClone("colz"); // merge histograms written by all threads
        cuts.DrawCut(layer);
        p++;
    }

    c11->cd(p);
    hEsilThetaLab.Merge()->DrawClone("colz");
    auto *gtheoKin{ActPhysics::Kinematics(TString::Format("%s(p,p)@65.2", beam.c_str()).Data()).GetKinematicLine3()};
    gtheoKin->SetLineColor(kMagenta);
    gtheoKin->Draw("l");

    auto *c12{new TCanvas{"c12", "Pipe1 PID L1"}};
    c12->DivideSquare(4);
    c12->cd(1);
    hl1.Merge()->DrawClone("colz");
    cuts.DrawCut("l1");
    c12->cd(2);
    hl1theta.Merge()->DrawClone("colz");
    c12->cd(3);
    hl1thetaCorr.Merge()->DrawClone("colz");
    auto *gtheo{ActPhysics::Kinematics(TString::Format("%s(p,p)@65.2", beam.c_str()).Data()).GetTheta3vs4Line()};
    gtheo->SetLineColor(46);
    gtheo->Draw("l");
    auto *gtheo2{ActPhysics::Kinematics(TString::Format("%s(p,4He)@65.2", beam.c_str()).Data()).GetTheta3vs4Line()};
    gtheo2->SetLineColor(3);
    gtheo2->Draw("l");
    auto l = new TLegend(0.6, 0.7, 0.9, 0.9);
    l->AddEntry(gtheo, gtheo->GetTitle());
    l->AddEntry(gtheo2, gtheo2->GetTitle());
    l->Draw();

    c12->cd(4);
    auto *c13{new TCanvas{"c13", "Pipe1 Temp"}};
    c13->cd();
    hl1.Merge()->DrawClone("colz");
    cuts.DrawCut("l1");
    cuts.DrawCut("l1p");
    cuts.DrawCut("l1bad");
}
