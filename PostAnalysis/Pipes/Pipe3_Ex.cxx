#ifndef Pipe3_Ex_cxx
#define Pipe3_Ex_cxx

#include "ActCutsManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe3_Ex(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    TString fIn;
    // if(light == "p")
    // fIn = TString::Format("./Outputs/tree_ESil_BSP_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    // else if(light == "4He")
    fIn = TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());

    // fIn = TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    // fIn = TString::Format("./Outputs/tree_pid_l1_p_diag2.root");


    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PID_Tree", fIn};

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
    std::string srimName {};
    if(light == "d")
        srimName = "2H";
    else if(light == "p")
        srimName = "1H";
    else if(light == "4He")
        srimName = "4He";
    srim->ReadTable(light,
                    TString::Format("../Simulation/SRIM/%s_H2-iC4H10_95-5_760mbar.txt", srimName.c_str()).Data());
    srim->ReadTable(beam, TString::Format("../Simulation/SRIM/%s_H2-iC4H10_95-5_760mbar.txt", beam.c_str()).Data());
    // srim->ReadTable("heavy", "../Simulation/SRIM/14O_H2-iC4H10_95-5_760mbar.txt");

    // Build energy at vertex
    auto dfVertex = df.Define("EVertex",
                              [&](const ActRoot::MergerData& d)
                              {
                                  double ret {};
                                  if(d.fLight.IsFilled())
                                      ret = srim->EvalInitialEnergy(light, d.fLight.fEs.front(), d.fLight.fTL);
                                  else // L1 trigger
                                      ret = srim->EvalEnergy(light, d.fLight.fTL);
                                  return ret;
                              },
                              {"MergerData"});

    // Init particles
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {target};
    ActPhysics::Particle pl {light};

    // Initial energy of beam at pad plane entrance
    double EBeamIni {3.84}; // MeV/u

    // // Filter on heavy particle hit in the telescope
    auto def {dfVertex};

    ActPhysics::Kinematics kin {pb, pt, pt, EBeamIni * pb.GetAMU()};
    auto beamThreshold {kin.GetT1Thresh()};
    // kin.Print();
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    std::vector<ActPhysics::Kinematics> vkins {def.GetNSlots()};
    for(auto& k : vkins)
        k = kin;

    double qval {3.923}; // qvalue of 17F+p
    auto lab_to_com {pt.GetMass() / (pb.GetMass() + pt.GetMass())};

    // If no Rec in title, quantity is found based on location of RPx + SRIM -> Uncertainty from initial beam energy,
    // RPx calculation and SRIM for beam. If Rec in title, quantity is reconstructed with missing mass method ->
    // uncertainty from EVertex (light particle energy, SRIM for proton), and angle.
    def = def.Define("EBeam",
                     [&](const ActRoot::MergerData& d)
                     {
                         auto ret {srim->Slow(beam, EBeamIni * pb.GetAMU(),
                                              d.fRP.X())}; // here RPx comes from merger data so it's in mm, not pads
                         if(ret <= 0)
                             ret = 1111111;
                         return ret;
                     },
                     {"MergerData"})
              .Define("ECM", [&](double EBeam) { return (EBeam * lab_to_com); }, {"EBeam"})
              .Define("ECN", [&](double Ecm) { return Ecm + qval; }, {"ECM"})
              .DefineSlot("RecEx",
                          [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                          {
                              if(EBeam < beamThreshold)
                                  return -1.;
                              else
                              {
                                  vkins[slot].SetBeamEnergy(EBeam);
                                  return vkins[slot].ReconstructExcitationEnergy(EVertex,
                                                                                 (d.fThetaLight) * TMath::DegToRad());
                              }
                          },
                          {"MergerData", "EVertex", "EBeam"})
              .DefineSlot("RecThetaCM",
                          [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
                          {
                              if(EBeam < beamThreshold)
                                  return -1.;
                              else
                              {
                                  vkins[slot].SetBeamEnergy(EBeam);
                                  return vkins[slot].ReconstructTheta3CMFromLab(EVertex,
                                                                                (d.fThetaLight) * TMath::DegToRad()) *
                                         TMath::RadToDeg();
                              }
                          },
                          {"MergerData", "EVertex", "EBeam"})
              .DefineSlot("RecEBeam",
                          [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex)
                          {
                              return vkins[slot].ReconstructBeamEnergyFromLabKinematics(EVertex, d.fThetaLight *
                                                                                                     TMath::DegToRad());
                          },
                          {"MergerData", "EVertex"})
              .Define("RecECM", [&](double rec_EBeam) { return lab_to_com * rec_EBeam; }, {"RecEBeam"})
              .Define("RecECN", [&](double rec_ECM) { return rec_ECM + qval; }, {"RecECM"});


    // // Book new histograms
    // auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    // auto hKinCM {def.Histo2D(HistConfig::KinCM, "ThetaCM", "EVertex")};

    // auto hEBeam {def.Histo1D(HistConfig::EBeam, "EBeam")};
    // auto hRecEBeam {def.Histo1D(HistConfig::EBeam, "RecEBeam")};
    // hRecEBeam->SetTitle("Reconstructed EBeam");

    // auto hExSil {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
    //                  .Histo1D(HistConfig::Ex, "Ex")};
    // hExSil->SetTitle("Ex with silicons");
    // auto hExL1 {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
    //                 .Histo1D(HistConfig::Ex, "Ex")};
    // hExL1->SetTitle("Ex with L1");

    // auto hTheta {def.Histo1D("fThetaLight")};

    // auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};

    // auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};
    // auto hRPx {def.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
    // auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "fThetaLight", "ThetaCM")};

    // auto hExESi {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
    //                  .Histo2D(HistConfig::EcmESi, "fLight.fEs", "ECM")};

    // // Ex dependences
    auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "RecThetaCM", "RecEx")};
    // auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    // auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};
    // auto hExSPz {def.Histo2D(HistConfig::ExZ, "fLight.fSP.fCoordinates.fZ", "Ex")};

    auto hEcnThetaCM {def.Histo2D(HistConfig::EcnThetaCM, "RecECN", "RecThetaCM")};
    // auto hEcn {def.Histo1D(HistConfig::ECN, "ECN")};
    // auto hEcm {def.Histo1D(HistConfig::ECM, "ECM")};
    // auto hECMRPX {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "ECM")};
    // // Beam energy against RP.X
    // auto hEBeamRPx {def.Histo2D({"hEBeamRPx", "Rec EBeam assuming elastic;RP.X [mm];EBeam", 300, 0, 260, 300, 0,
    // 300},
    //                             "fRP.fCoordinates.fX", "RecEBeam")};

    // Make histograms

    // Ex Si
    auto hExESi {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                     .Histo2D(HistConfig::ExESi, "fLight.fEs", "RecEx")};
    auto hRecExSil {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                        .Histo1D(HistConfig::Ex, "RecEx")};
    hRecExSil->SetTitle("Reconstructed Ex with silicons");
    auto hRecExSilRPx {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                           .Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "RecEx")};
    hRecExSilRPx->GetYaxis()->SetRangeUser(-5., 5.);
    hRecExSilRPx->SetTitle("Rec Ex with Si vs RPx");
    auto hRecExL1 {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                       .Histo1D(HistConfig::Ex, "RecEx")};
    hRecExL1->SetTitle("Reconstructed Ex with L1");
    auto hRecExL1RPx {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                          .Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "RecEx")};
    hRecExL1RPx->GetYaxis()->SetRangeUser(-5., 5.);
    hRecExL1RPx->SetTitle("Rec Ex with L1 vs RPx");


    // Ebeam
    auto hEbeam {def.Histo1D(HistConfig::EBeam, "EBeam")};
    auto hRecEbeam {def.Histo1D(HistConfig::EBeam, "RecEBeam")};
    auto hEBeamCompare {def.Histo2D(HistConfig::EBeamCompare, "EBeam", "RecEBeam")};
    auto hEBeamRPx {def.Histo2D(HistConfig::EBeamRPx, "fRP.fCoordinates.fX", "EBeam")};
    hEBeamRPx->SetTitle("EBeam from srim vs RPx");
    auto hRecEBeamRPx {def.Histo2D(HistConfig::EBeamRPx, "fRP.fCoordinates.fX", "RecEBeam")};
    hRecEBeamRPx->SetTitle("Reconstructed EBeam vs RPx");

    // EVertex
    auto hEVertexRPx {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "EVertex")};
    hEVertexRPx->SetTitle("EVertex vs RPx; RPx [mm]; EVertex [MeV]");
    hEVertexRPx->GetYaxis()->SetRangeUser(0., 15.);
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};
    auto hKinCM {def.Histo2D(HistConfig::KinCM, "RecThetaCM", "EVertex")};
    auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "RecThetaCM", "fThetaLight")};
    auto hRPxThetaCM {def.Histo2D(HistConfig::RPxThetaCM, "fRP.fCoordinates.fX", "RecThetaCM")};


    // Ecm
    auto hEcnSi {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                     .Histo1D(HistConfig::ECN, "ECN")};

    auto hRecEcnSi {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                        .Histo1D(HistConfig::ECN, "RecECN")};
    auto hRecEcnSiFront {def.Filter([](ActRoot::MergerData& m)
                                    { return m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"; },
                                    {"MergerData"})
                             .Histo1D(HistConfig::ECN, "RecECN")};
    hRecEcnSiFront->SetTitle("Front");
    auto hRecEcnSiLat {
        def.Filter(
               [](ActRoot::MergerData& m)
               { return m.fLight.GetNLayers() == 1 && (m.fLight.GetLayer(0) == "r0" || m.fLight.GetLayer(0) == "l0"); },
               {"MergerData"})
            .Histo1D(HistConfig::ECN, "RecECN")};
    hRecEcnSiLat->SetTitle("Lat");

    auto hEcnL1 {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                     .Histo1D(HistConfig::ECN, "ECN")};
    auto hRecEcnL1 {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                        .Histo1D(HistConfig::ECN, "RecECN")};
    hRecEcnSi->SetTitle("Ecn from Si");
    hRecEcnL1->SetTitle("Ecn from L1");

    // Draw
    auto* c31 {new TCanvas("c31", "Pipe 3 Canvas 1: Ex", 1000, 800)};
    c31->DivideSquare(6);
    c31->cd(1);
    // gPad->SetLogy();
    hRecExSil->DrawClone();
    c31->cd(2);
    hRecExSilRPx->DrawClone("colz");
    c31->cd(3);
    hExESi->DrawClone("colz");
    c31->cd(4);
    gPad->SetLogy();
    hRecExL1->DrawClone();
    c31->cd(5);
    hRecExL1RPx->DrawClone("colz");
    c31->cd(6);
    hEcnThetaCM->DrawClone("colz");

    TLine* hgs = new TLine(0, 0, 0, 50000);
    TLine* h1ex = new TLine(0.495, 0, 0.495, 50000);
    hgs->SetLineColor(2);
    h1ex->SetLineColor(2);
    hgs->SetLineStyle(2);
    h1ex->SetLineStyle(2);
    c31->cd(1);
    hgs->Draw("same");
    h1ex->Draw("same");
    c31->cd(4);
    hgs->Draw("same");
    h1ex->Draw("same");

    TLine* vgs = new TLine(0, 0, 50000, 0);
    TLine* v1ex = new TLine(0, 0.495, 50000, 0.495);
    TLine* v2ex = new TLine(0, 3.104, 50000, 3.104);
    vgs->SetLineColor(2);
    v1ex->SetLineColor(2);
    v2ex->SetLineColor(2);
    vgs->SetLineStyle(2);
    v1ex->SetLineStyle(2);
    v2ex->SetLineStyle(2);
    c31->cd(3);
    vgs->Draw("same");
    v1ex->Draw("same");
    v2ex->Draw("same");

    auto* c32 {new TCanvas("c32", "Pipe 3 Canvas 2: EBeam", 800, 800)};
    c32->DivideSquare(4);
    c32->cd(1);
    hEBeamRPx->DrawClone("colz");
    c32->cd(2);
    hRecEBeamRPx->DrawClone("colz");
    c32->cd(3);
    hRecEbeam->DrawClone();
    hEbeam->SetLineColor(2);
    hEbeam->DrawClone("same");
    auto t1 = new TLatex(40, 700, "reconstructed");
    t1->SetTextColor(1);
    t1->Draw("same");
    auto t2 = new TLatex(40, 650, "from srim");
    t2->SetTextColor(2);
    t2->Draw("same");
    c32->cd(4);
    hEBeamCompare->DrawClone("colz");

    auto* c33 {new TCanvas("c33", "Pipe 3 Canvas 3: EVertex", 800, 600)};
    c33->DivideSquare(4);
    c33->cd(1);
    hEVertexRPx->DrawClone("colz");
    c33->cd(2);
    hRPxThetaCM->DrawClone("colz");
    c33->cd(3);
    hKinCM->DrawClone("colz");
    c33->cd(4);
    hThetaCMLab->DrawClone("colz");
    auto* gThetaLabvsCM {kin.GetThetaLabvsThetaCMLine()};
    gThetaLabvsCM->Draw("same");


    auto* c34 {new TCanvas("c34", "Pipe 3 Canvas 4: ECN", 1500, 600)};
    c34->DivideSquare(2);
    c34->cd(1);
    // hRecEcnSi->DrawClone();
    // hRecEcnSiLat->SetLineColor(2);
    // hRecEcnSiLat->DrawClone("same");
    // hRecEcnSiFront->SetLineColor(8);
    // hRecEcnSiFront->DrawClone("same");
    // c34->BuildLegend();

    auto htot = (TH1*)hRecEcnSi->DrawClone();
    hRecEcnSiLat->SetLineColor(2);
    auto hlat = (TH1*)hRecEcnSiLat->DrawClone("same");
    hRecEcnSiFront->SetLineColor(kMagenta);
    auto hfront = (TH1*)hRecEcnSiFront->DrawClone("same");
    auto leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(htot, "All", "l");
    leg->AddEntry(hlat, "Lateral", "l");
    leg->AddEntry(hfront, "Front", "l");
    leg->Draw();

    // hEcnSi->SetLineColor(2);
    // hEcnSi->DrawClone("same");
    // auto t41 = new TLatex(6, 700, "reconstructed");
    // t41->SetTextColor(1);
    // t41->Draw("same");
    // auto t42 = new TLatex(6, 620, "from srim");
    // t42->SetTextColor(2);
    // t42->Draw("same");
    c34->cd(2);
    hRecEcnL1->DrawClone();
    // hEcnL1->SetLineColor(2);
    // hEcnL1->DrawClone("same");

    //   States to draw
    std::vector<double> states18Ne {{4517., 4524., 4590., 5100., 5142., 5461., 6135., 6301., 6352., 6850., 7050., 7120.,
                                     7352., 7601., 7717., 7924.}};
    for(auto& s : states18Ne)
    {
        TLine* st = new TLine(s * 1e-3, 0, s * 1e-3, 900);
        st->SetLineColorAlpha(16,0.2);
        st->SetLineStyle(2);
        c34->cd(1);
        st->Draw("same");
        c34->cd(2);
        // st->Draw("same");
    }

    // Save!
    // auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s_phi_lt2.root", beam.c_str(), target.c_str(),
    // light.c_str())};
    auto outfile {
        TString::Format("./Outputs/tree_ex_%s_%s_%s_%.2f.root", beam.c_str(), target.c_str(), light.c_str(), EBeamIni)};
    def.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';


    // ActRoot::CutsManager<std::string> cuts;
    // cuts.ReadCut("cut", TString::Format("./Cuts/eVertexRpx_side_wall2.root").Data());
    // TH1D* hangles {new TH1D("hangles", "hangles", 100, 0, 180)};
    // auto listOfCuts {cuts.GetListOfKeys()};
    // if(listOfCuts.size())
    // {
    //     // Apply cut and save in file
    //     auto gated {def.Filter(
    //         [&](double evertex, ActRoot::MergerData& m)
    //         {
    //             if(cuts.GetCut("cut"))
    //             {
    //                 return cuts.IsInside("cut", m.fRP.X(), evertex);
    //             }
    //             else
    //                 return false;
    //         },
    //         {"EVertex", "MergerData"})};
    //     auto name {TString::Format("./Outputs/tree_Ex_side_walls2.root")};
    //     std::cout << "Saving in file : " << name << '\n';
    //     gated.Snapshot("Final_Tree", name.Data());

    //     // Draw the cut
    //     c33->cd(1);
    //     cuts.DrawCut("cut");
    //     auto* c30 {new TCanvas("c35", "Pipe 3 Canvas 0: testing", 800, 600)};
    //     c30->cd();
    //     gated.Foreach([&](ActRoot::MergerData& m) { hangles->Fill(m.fThetaLight); }, {"MergerData"});
    //     hangles->Draw();
    // }
}
#endif
