#ifndef Pipe2_Ex_cxx
#define Pipe2_Ex_cxx

#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../HistConfig.h"

void Pipe2_Ex(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PID_Tree", filename};

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
              .Define("ECN", [&](double EBeam) { return (EBeam * lab_to_com + qval); }, {"EBeam"})
              .Define("ECM", [&](double EBeam) { return (EBeam * lab_to_com); }, {"EBeam"})
              .DefineSlot("Ex",
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
              .DefineSlot("ThetaCM",
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
              .Define("Rec_ECM", [&](double rec_EBeam) { return lab_to_com * rec_EBeam; }, {"RecEBeam"})
              .Define("BSP",
              [&](const ActRoot::MergerData& m)
              {
                  double bsp = m.fRP.X() + m.fHeavy.fTL;
                  return bsp;
              },
              {"MergerData"});


    // Book new histograms
    auto hKin {def.Histo2D(HistConfig::KinEl, "fThetaLight", "EVertex")};

    auto hKinCM {def.Histo2D(HistConfig::KinCM, "ThetaCM", "EVertex")};

    auto hEBeam {def.Histo1D(HistConfig::EBeam, "EBeam")};
    auto hRecEBeam {def.Histo1D(HistConfig::EBeam, "RecEBeam")};
    hRecEBeam->SetTitle("Reconstructed EBeam");
    auto hExSil {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                     .Histo1D(HistConfig::Ex, "Ex")};
    hExSil->SetTitle("Ex with silicons");
    auto hExL1 {def.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == false; }, {"MergerData"})
                    .Histo1D(HistConfig::Ex, "Ex")};
    hExL1->SetTitle("Ex with L1");

    auto hTheta {def.Histo1D("fThetaLight")};

    auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};

    auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};
    auto hRPx {def.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
    auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "fThetaLight", "ThetaCM")};

    auto hExESi {def.Filter([](ActRoot::MergerData& m) { return m.fLight.GetNLayers() == 1 ; }, {"MergerData"}).Histo2D(HistConfig::EcmESi, "fLight.fEs","ECM")};

    // Ex dependences
    auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "ThetaCM", "Ex")};
    auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};
    auto hExSPz {def.Histo2D(HistConfig::ExZ, "fLight.fSP.fCoordinates.fZ", "Ex")};

    auto hEcnThetaCM {def.Histo2D(HistConfig::EcnThetaCM, "ECN", "ThetaCM")};
    auto hEcn {def.Histo1D(HistConfig::ECN, "ECN")};
    auto hEcm {def.Histo1D(HistConfig::ECM, "ECM")};
    auto hECMRPX {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "ECM")};
    // Beam energy against RP.X
    auto hEBeamRPx {def.Histo2D({"hEBeamRPx", "Rec EBeam assuming elastic;RP.X [mm];EBeam", 300, 0, 260, 300, 0, 300},
                                "fRP.fCoordinates.fX", "RecEBeam")};

    // Save!
    auto outfile {TString::Format("./Outputs/tree_ex_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    def.Snapshot("Final_Tree", outfile);
    std::cout << "Saving Final_Tree in " << outfile << '\n';
    // Save Ex histos on file
    auto file {std::make_shared<TFile>(outfile.Data(), "update")};
    hExSil->Write("hExSil");
    hExL1->Write("hExL1");
    file->Close();


    auto* c22 {new TCanvas("c22", "Pipe2 canvas 2")};
    c22->DivideSquare(4);
    c22->cd(1);
    hRP->DrawClone("colz");
    c22->cd(2);
    hExESi->DrawClone("colz");
    c22->cd(3);
    hEBeam->DrawClone();
    c22->cd(4);
    hRecEBeam->DrawClone();

    auto* c21 {new TCanvas("c21", "Pipe2 canvas 1")};
    c21->DivideSquare(6);
    c21->cd(1);
    hKin->DrawClone("colz");
    auto* theo {kin.GetKinematicLine3()};
    theo->Draw("same");
    c21->cd(2);
    hExSil->DrawClone();
    c21->cd(3);
    hExL1->DrawClone();
    // hKinCM->DrawClone("colz");
    c21->cd(4);
    // hExL1->DrawClone();
    // hExThetaLab->DrawClone("colz");
    hExSPz->DrawClone("colz");
    c21->cd(5);
    hExThetaCM->DrawClone("colz");
    c21->cd(6);
    hExRP->DrawClone("colz");

    auto* c23 {new TCanvas {"c23", "Pipe2 canvas 3"}};
    c23->DivideSquare(6);
    c23->cd(1);
    hThetaCMLab->DrawClone("colz");
    c23->cd(2);
    hEBeamRPx->DrawClone("colz");
    c23->cd(3);
    hEcnThetaCM->DrawClone("colz");
    c23->cd(4);
    hEcn->DrawClone();
    c23->cd(5);
    hEcm->DrawClone();
    c23->cd(6);
    hECMRPX->DrawClone("colz");
}
#endif
