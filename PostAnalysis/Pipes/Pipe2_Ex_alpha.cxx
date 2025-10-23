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

void Pipe2_Ex_alpha(const std::string& beam, const std::string& target, const std::string& light)
{
    // Read data
    auto filename {TString::Format("./Outputs/tree_pid_%s_%s_%s.root", beam.c_str(), target.c_str(), light.c_str())};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {"PID_Tree", filename};

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
    std::string srimName {};
    srim->ReadTable(light, TString::Format("../Simulation/SRIM/%s_H2-iC4H10_95-5_760mbar.txt", light.c_str()).Data());
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

    // ActPhysics::Kinematics kin {pb, pt, pt, EBeamIni * pb.GetAMU()};
    // auto beamThreshold {kin.GetT1Thresh()};
    // kin.Print();
    // Vector of kinematics as one object is needed per
    // processing slot (since we are changing EBeam in each entry)
    // std::vector<ActPhysics::Kinematics> vkins {def.GetNSlots()};
    // for(auto& k : vkins)
    //     k = kin;

    double qval {3.923}; // qvalue of 17F+p
    auto lab_to_cm {pt.GetMass() / (pb.GetMass() + pt.GetMass())};

    def = def.Define("EBeam",
                     [&](const ActRoot::MergerData& m)
                     {
                         auto ret {srim->Slow(beam, EBeamIni * pb.GetAMU(),
                                              m.fRP.X())}; // here RPx comes from merger data so it's in mm, not pads
                         if(ret <= 0)
                             ret = 1111111;
                         return ret;
                     },
                     {"MergerData"})
              .Define("ECN", [&](double EBeam) { return (EBeam * lab_to_cm + qval); }, {"EBeam"})
              .Define("ECM", [&](double EBeam) { return (EBeam * lab_to_cm); }, {"EBeam"})
              .Define("lightE",
                      [&](const ActRoot::MergerData& m)
                      {
                          auto ret {srim->EvalEnergy(light, m.fLight.fTL)};
                          return ret;
                      },
                      {"MergerData"})
              .Define("heavyE",
                      [&](const ActRoot::MergerData& m)
                      {
                          auto ret {srim->EvalEnergy("heavy", m.fHeavy.fTL)};
                          return ret;
                      },
                      {"MergerData"})
                .Define("RecECN",[&](double l, double h){return l+h;},{"lightE","heavyE"});
    //   .DefineSlot("Ex",
    //               [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
    //               {
    //                   if(EBeam < beamThreshold)
    //                       return -1.;
    //                   else
    //                   {
    //                       vkins[slot].SetBeamEnergy(EBeam);
    //                       return vkins[slot].ReconstructExcitationEnergy(EVertex,
    //                                                                      (d.fThetaLight) * TMath::DegToRad());
    //                   }
    //               },
    //               {"MergerData", "EVertex", "EBeam"})
    //   .DefineSlot("ThetaCM",
    //               [&](unsigned int slot, const ActRoot::MergerData& d, double EVertex, double EBeam)
    //               {
    //                   if(EBeam < beamThreshold)
    //                       return -1.;
    //                   else
    //                   {
    //                       vkins[slot].SetBeamEnergy(EBeam);
    //                       return vkins[slot].ReconstructTheta3CMFromLab(EVertex,
    //                                                                     (d.fThetaLight) * TMath::DegToRad()) *
    //                              TMath::RadToDeg();
    //                   }
    //               },
    //               {"MergerData", "EVertex", "EBeam"});


    // Book new histograms

    auto hEBeam {def.Histo1D(HistConfig::EBeam, "EBeam")};
    // auto hEx {def.Histo1D(HistConfig::Ex, "Ex")};

    auto hTheta {def.Histo1D("fThetaLight")};

    auto hThetaBeam {def.Histo2D(HistConfig::ThetaBeam, "fRP.fCoordinates.fX", "fThetaBeam")};

    auto hRP {def.Histo2D(HistConfig::RP, "fRP.fCoordinates.fX", "fRP.fCoordinates.fY")};
    auto hRPx {def.Histo1D(HistConfig::RPx, "fRP.fCoordinates.fX")};
    // auto hThetaCMLab {def.Histo2D(HistConfig::ThetaCMLab, "fThetaLight", "ThetaCM")};

    // Ex dependences
    // auto hExThetaCM {def.Histo2D(HistConfig::ExThetaCM, "ThetaCM", "Ex")};
    // auto hExThetaLab {def.Histo2D(HistConfig::ExThetaLab, "fThetaLight", "Ex")};
    // auto hExRP {def.Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "Ex")};
    // auto hExSPz {def.Histo2D(HistConfig::ExZ, "fLight.fSP.fCoordinates.fZ", "Ex")};

    // auto hEcnThetaCM {def.Histo2D(HistConfig::EcnThetaCM, "ECN", "ThetaCM")};
    auto hEcn {def.Histo1D(HistConfig::ECN, "ECN")};
    auto hEcm {def.Histo1D(HistConfig::ECM, "ECM")};
    auto hECMRPX {def.Histo2D(HistConfig::RPxECM, "fRP.fCoordinates.fX", "ECM")};
    auto hEalpha {def.Histo1D(HistConfig::ELight, "lightE")};
    auto hEheavy {def.Histo1D(HistConfig::EHeavy, "heavyE")};
    auto hEcnRec {def.Histo1D(HistConfig::ECN, "RecECN")};

    hEheavy->SetTitle("14O energy");

    auto* c22 {new TCanvas("c22", "Pipe2 canvas 2")};   
    c22->DivideSquare(6);
    c22->cd(1);
    hRP->DrawClone("colz");
    c22->cd(2);
    hRPx->DrawClone();
    c22->cd(3);
    hEBeam->DrawClone();
    c22->cd(4);
    hEcm->DrawClone();
    c22->cd(5);
    hEcn->DrawClone();
    hEcnRec->SetLineColor(2);
    // hEcnRec->DrawClone("same");
    c22->cd(6);
    hECMRPX->DrawClone("colz");

    auto* c21 {new TCanvas("c21", "Pipe2 canvas 1")};
    c21->DivideSquare(4);
    c21->cd(1);
    hEalpha->DrawClone();
    c21->cd(2);
    hEheavy->DrawClone();
}
#endif
