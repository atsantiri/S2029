#include "ActKinematics.h"
#include "ActParticle.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TLine.h"

#include "../../PostAnalysis/HistConfig.h"
#include "../Classes/DoubleXS.cxx"
#include "../Classes/DoubleXS.h"

void Run()
{
    ROOT::EnableImplicitMT();
    double EBeamIni {3.84};
    double pressure {760};
    std::string beam {"17F"};
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {"1H"};

    // Pipe 3 output
    ROOT::RDataFrame df {"Final_Tree", TString::Format("../../PostAnalysis/Outputs/tree_ex_%s_p_p_%.02f_sil.root",
                                                       beam.c_str(), EBeamIni)};
    auto hEx {
        df.Histo1D({"hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (4. - (-2.)) / 200 * 1e3),
                    200, -2, 4},
                   "RecEx")};
    auto h2DECN {df.Histo2D(HistConfig::EcnThetaCM, "RecECN", "RecThetaCM")};
    h2DECN->SetTitle("Experiment");
    auto h2DECM {df.Histo2D(HistConfig::Eff2D, "RecECM", "RecThetaCM")};
    h2DECM->SetTitle("Experiment");

    // Simulation output
    auto simuFile {std::make_unique<TFile>(
        TString::Format("../../Simulation/Outputs/Simu_%s_p_p_%.0fmbar.root", beam.c_str(), pressure))};
    auto heff {simuFile->Get<TH2D>("hEff2D")};
    heff->SetTitle("Simu eff");
    heff->SetDirectory(nullptr);
    auto heffCN {simuFile->Get<TH2D>("hEffECN2D")};
    heffCN->SetTitle("Simu eff E_{^{18}Ne}");
    heffCN->SetDirectory(nullptr);
    simuFile->Close();

    // Read srim
    auto srim {new ActPhysics::SRIM};
    srim->ReadTable("beam", TString::Format("../../Simulation/SRIM/%s_H2-iC4H10_95-5_%.0fmbar.txt", beam.c_str(),pressure).Data());

    // Kinematics
    // ActPhysics::Kinematics kin {pb, pt, pt, EBeamIni * pb.GetAMU()};
    auto kin {new ActPhysics::Kinematics {TString::Format("%s(p,p)@%.0f",beam.c_str(),EBeamIni*pb.GetAMU()).Data()}};

    // Number of beam particles from Pipe_Beam
    double Nb {176922 * 300}; // counter of gatconf == 64 * div factor

    // Target density: tbd from LISE++
    double rho {1.};


    // Draw
    auto* c0 {new TCanvas("c0", "Canvas for inspection 0", 1500, 1000)};
    c0->DivideSquare(6);
    c0->cd(1);
    h2DECN->DrawClone("colz");
    c0->cd(2);
    heffCN->DrawClone("colz");
    c0->cd(3);
    // hnorm->DrawClone("colz");
    c0->cd(4);
    h2DECM->DrawClone("colz");
    c0->cd(5);
    heff->DrawClone("colz");
    c0->cd(6);
    // hnormCN->DrawClone("colz");

    DoubleXS xs {h2DECM.GetPtr(), heff, srim, Nb, rho, kin};
    xs.Draw();
}