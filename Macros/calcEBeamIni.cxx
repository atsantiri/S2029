///////////////////////////////////////////////////////////////////////////////////////////
// script to calculate average energy of 17F beam in entrance of ACTAR from the          //
// track lengths                                                                         //
///////////////////////////////////////////////////////////////////////////////////////////

#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TROOT.h"


void calcEBeamIni()
{
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EMerge};
    // dataman.SetRuns(49, 65);
    dataman.SetRuns(49, 49);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df {*chain};


    auto gated =
        df.Filter([](ActRoot::TPCData& tpc, ActRoot::ModularData& mod)
                  { return (mod.Get("GATCONF") == 64) && (tpc.fClusters.size() == 1); }, {"TPCData", "ModularData"})
            .Define("fLxy",
                    [&](ActRoot::TPCData& tpc)
                    {
                        auto cluster {tpc.fClusters[0]};
                        auto line {cluster.GetRefToLine()};
                        auto dir {line.GetDirection()};
                        cluster.SortAlongDir(dir);
                        auto lastVoxel {cluster.GetRefToVoxels().back()};
                        auto firstVoxel {cluster.GetRefToVoxels().front()};
                        auto lastPoint {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
                        auto firstPoint {line.ProjectionPointOnLine(firstVoxel.GetPosition())};
                        double lxy = TMath::Sqrt(TMath::Power(lastPoint.X() - firstPoint.X(), 2) +
                                                 TMath::Power(lastPoint.Y() - firstPoint.Y(), 2));
                        return lxy * 2; // Conversion factor from pads to mm
                    },
                    {"TPCData"});

    ROOT::TThreadedObject<TH2D> h2d {"hPad", "Pad plane ;X [pad];Y [pad]", 128, 0, 128, 128, 0, 128};
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
            }
        },
        {"TPCData"});

    auto* c1 {new TCanvas {"c1", "c1", 1400, 1400}};
    c1->DivideSquare(4);
    c1->cd(1);
    h2d->DrawClone("colz");

    ROOT::TThreadedObject<TH1D> hEini {"hEini", "beam E initial; ELab [MeV/u];Counts", 100, 2.5, 4.5};
    ROOT::TThreadedObject<TH1D> hres {"hres", "Resonance location; x [mm];Counts", 100, 0, 256};

    auto* srim {new ActPhysics::SRIM()};
    ActPhysics::Particle beam {"17F"};
    ActPhysics::Particle target {"p"};
    srim->ReadTable("beamInGas", "../Simulation/SRIM/17F_H2-iC4H10_95-5_760mbar.txt");

    auto qval {3.9231};
    auto resonance {6.150};
    auto LStoCMS {target.GetAMU()/(target.GetAMU()+beam.GetAMU())};
    auto Eres = (resonance-qval)/LStoCMS;
    std::cout<<"resonance energy in LS: "<<Eres/beam.GetAMU()<<" MeV/u"<<std::endl;

    gated.Foreach(
        [&](double lxy)
        {
            auto ene {srim->EvalInitialEnergy("beamInGas", 0, lxy)};
            hEini->Fill(ene / beam.GetAMU());

            auto dist {srim->TravelledDistance("beamInGas", ene, Eres)};
            // std::cout<<ene<<" "<<Eres<<" "<<dist<<std::endl;
            hres->Fill(dist);
        },
        {"fLxy"});

    c1->cd(2);
    TF1* fE = new TF1("fE", "gaus(0)", 3.5, 4);
    fE->SetParameters(2400, 3.8, 0.06);
    hEini->Fit(fE, "Q");
    auto meanEne {fE->GetParameter(1)};
    std::cout << "Initial Beam Energy: " << std::fixed << std::setprecision(3) << meanEne << " +- "
              << 2.355 * fE->GetParameter(2) << " (" << 2.355 * fE->GetParameter(2) / meanEne * 100
              << " %) or in CMS ECM = "<< meanEne*beam.GetAMU()*LStoCMS << std::endl;
    hEini->DrawClone();

    c1->cd(3);
    TF1* fdist = new TF1("fE", "gaus(0)", 120, 160);
    fdist->SetParameters(2400, 150, 0.06);
    hres->Fit(fdist, "Q");
    auto meanRP {fdist->GetParameter(1)};
    std::cout << "Resonance will be around: " << std::fixed << std::setprecision(3) << fdist->GetParameter(1) << " +- "
              << 2.355 * fdist->GetParameter(2) << std::endl;
    hres->DrawClone();

}