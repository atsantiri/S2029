#ifndef Pipe3_Ex_cxx
#define Pipe3_Ex_cxx

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
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
    if(light == "p")
        fIn = TString::Format("./Outputs/tree_ESil_BSP_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());
    else if(light == "4He")
        fIn = TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str());

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
    srim->ReadTable("heavy", "../Simulation/SRIM/14O_H2-iC4H10_95-5_760mbar.txt");

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
    kin.Print();
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
              .Define("RecECM", [&](double rec_EBeam) { return lab_to_com * rec_EBeam; }, {"RecEBeam"});


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
    // auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "ThetaCM", "Ex")};
    // auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    // auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};
    // auto hExSPz {def.Histo2D(HistConfig::ExZ, "fLight.fSP.fCoordinates.fZ", "Ex")};

    // auto hEcnThetaCM {def.Histo2D(HistConfig::EcnThetaCM, "ECN", "ThetaCM")};
    // auto hEcn {def.Histo1D(HistConfig::ECN, "ECN")};
    // auto hEcm {def.Histo1D(HistConfig::ECM, "ECM")};
    // auto hECMRPX {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "ECM")};
    // // Beam energy against RP.X
    // auto hEBeamRPx {def.Histo2D({"hEBeamRPx", "Rec EBeam assuming elastic;RP.X [mm];EBeam", 300, 0, 260, 300, 0,
    // 300},
    //                             "fRP.fCoordinates.fX", "RecEBeam")};

    // Make histograms
    auto hExESi {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                     .Histo2D(HistConfig::ExESi, "fLight.fEs", "RecEx")};
    auto hRecExSil {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "RecEx")};
    hRecExSil->SetTitle("Reconstructed Ex with silicons");


    auto hEbeam {def.Histo1D(HistConfig::EBeam, "EBeam")};
    auto hRecEbeam {def.Histo1D(HistConfig::EBeam, "RecEBeam")};
    auto hEBeamCompare {def.Histo2D(HistConfig::EBeamCompare, "EBeam", "RecEBeam")};
    auto hEBeamRPx {def.Histo2D(HistConfig::EBeamRPx, "fRP.fCoordinates.fX", "EBeam")};
    hEBeamRPx->SetTitle("EBeam from srim vs RPx");
    auto hRecEBeamRPx {def.Histo2D(HistConfig::EBeamRPx, "fRP.fCoordinates.fX", "RecEBeam")};
    hRecEBeamRPx->SetTitle("Reconstructed EBeam vs RPx");

    // Draw
    auto* c1 {new TCanvas("c1", "Ex", 800, 800)};
    c1->DivideSquare(3);
    c1->cd(1);
    hRecExSil->DrawClone();

    c1->cd(2);
    hExESi->DrawClone("colz");


    auto* c2 {new TCanvas("c2", "EBeam", 800, 800)};
    c2->DivideSquare(4);
    c2->cd(1);
    hEBeamRPx->DrawClone("colz");
    c2->cd(2);
    hRecEBeamRPx->DrawClone("colz");
    c2->cd(3);
    hRecEbeam->DrawClone();
    hEbeam->SetLineColor(2);
    hEbeam->DrawClone("same");
    auto t1 = new TLatex(40, 700, "reconstructed");
    t1->SetTextColor(1);
    t1->Draw("same");
    auto t2 = new TLatex(40, 650, "from srim");
    t2->SetTextColor(2);
    t2->Draw("same");
    c2->cd(4);
    hEBeamCompare->DrawClone("colz");

    // Save!
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    def.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';
}
#endif
