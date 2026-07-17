#include "ActKinematics.h"
#include "ActParticle.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TLine.h"

#include "../HistConfig.h"

void normalizeEx()
{
    double EBeamIni {3.84};
    std::string beam {"17F"};
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {"1H"};

    // Pipe 3 output
    ROOT::RDataFrame df {"Final_Tree",
                         TString::Format("../Outputs/tree_ex_%s_p_p_%.02f_sil.root", beam.c_str(), EBeamIni)};

    // Simulation output
    auto simuFile {
        std::make_unique<TFile>(TString::Format("../../Simulation/Outputs/Simu_%s_p_p_760mbar.root", beam.c_str()))};
    auto heff {simuFile->Get<TH2D>("hEff2D")};
    heff->SetTitle("Simu eff");
    heff->SetDirectory(nullptr);
    auto heffCN {simuFile->Get<TH2D>("hEffECN2D")};
    heffCN->SetTitle("Simu eff E_{^{18}Ne}");
    heffCN->SetDirectory(nullptr);
    simuFile->Close();

    // Kinematics
    ActPhysics::Kinematics kin {pb, pt, pt, EBeamIni * pb.GetAMU()};

    // Number of beam particles from Pipe_Beam
    double Nb {176922 * 300}; // counter of gatconf == 64 * div factor

    // Target density: tbd from LISE++
    double rho {1.};

    // Book histograms
    auto hEx {
        df.Histo1D({"hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (4. - (-2.)) / 200 * 1e3),
                    200, -2, 4},
                   "RecEx")};
    auto h2DECN {df.Histo2D(HistConfig::EcnThetaCM, "RecECN", "RecThetaCM")};
    h2DECN->SetTitle("Experiment");
    auto h2DECM {df.Histo2D(HistConfig::Eff2D, "RecECM", "RecThetaCM")};
    h2DECM->SetTitle("Experiment");

    // Normalize
    auto hnorm = (TH2D*)h2DECM->Clone("hnorm");
    hnorm->SetTitle("Efficiency corrected");
    hnorm->Sumw2();
    heff->Sumw2();
    hnorm->Divide(heff);

    auto hnormCN = (TH2D*)h2DECN->Clone("hnormCN");
    hnormCN->SetTitle("Efficiency Corrected CN");
    hnormCN->Sumw2();
    heffCN->Sumw2();
    hnormCN->Divide(heffCN);

    // Draw
    auto* c0 {new TCanvas("c0", "Canvas for inspection 0", 1500, 1000)};
    c0->DivideSquare(6);
    c0->cd(1);
    h2DECN->DrawClone("colz");
    c0->cd(2);
    heffCN->DrawClone("colz");
    c0->cd(3);
    hnorm->DrawClone("colz");
    c0->cd(4);
    h2DECM->DrawClone("colz");
    c0->cd(5);
    heff->DrawClone("colz");
    c0->cd(6);
    hnormCN->DrawClone("colz");


    TH1D* hECN = hnormCN->ProjectionX("hCN");
    auto* c1 {new TCanvas("c1", "Canvas for inspection 1", 900, 600)};
    hECN->DrawClone("hist");
}