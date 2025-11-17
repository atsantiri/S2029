#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"


void beamBraggPeak()
{

    bool compare_with_simu {true};

    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    dataman.SetRuns(49, 49);
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    auto gated {
        df.Filter([](ActRoot::TPCData& tpc, ActRoot::ModularData& mod)
                  { return (mod.Get("GATCONF") == 64) && (tpc.fClusters.size() == 1); }, {"TPCData", "ModularData"})
            .Define(
                "lastPoint",
                [](ActRoot::TPCData& tpc)
                {
                    auto cl {tpc.fClusters[0]};              // Get the cluster
                    auto l {cl.GetRefToLine()};              // Get the line that was fitted on the cluster
                    auto dir {l.GetDirection()};             // Get its direction
                    cl.SortAlongDir(dir);                    // sort along the line to find first/last voxel
                    auto voxel {cl.GetRefToVoxels().back()}; // get voxel of interest - here the last
                    auto projection {l.ProjectionPointOnLine(
                        voxel.GetPosition())}; // Get the projection of the position of the first/last voxel on the line
                    return projection;
                },
                {"TPCData"})
            .Define(
                "firstPoint",
                [](ActRoot::TPCData& tpc)
                {
                    auto cl {tpc.fClusters[0]};               // Get the cluster
                    auto l {cl.GetRefToLine()};               // Get the line that was fitted on the cluster
                    auto dir {l.GetDirection()};              // Get its direction
                    cl.SortAlongDir(dir);                     // sort along the line to find first/last voxel
                    auto voxel {cl.GetRefToVoxels().front()}; // get voxel of interest - here the first
                    auto projection {l.ProjectionPointOnLine(
                        voxel.GetPosition())}; // Get the projection of the position of the first/last voxel on the line
                    return projection;
                },
                {"TPCData"})
            .Define("lastPY", "lastPoint.Y()")
            .Define("firstPY", "firstPoint.Y()")
            .Define("dY", "lastPY-firstPY")
            .Filter([](float fPY, float lPY, float dY)
                    { return (fPY >= 59 && fPY <= 64) && (lPY >= 59 && lPY <= 64) && (dY >= -1 && dY <= 1); },
                    {"firstPY", "lastPY", "dY"}) // constrain to forward angles
    };

    auto h2d = new TH2D {"hPad", "Pad plane;X [pad];Y [pad]", 128, 0, 128, 128, 0, 128};
    auto hBragg = new TH1D {"hBragg", "; X [pad]; Q per pad", 128, 0, 128};
    gated.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            auto& cluster {tpc.fClusters.front()};
            for(const auto& v : cluster.GetVoxels())
            {
                auto posX {v.GetPosition().X()};
                auto posY {v.GetPosition().Y()};
                auto q {v.GetCharge()};
                hBragg->Fill(posX, q);
                h2d->Fill(posX, posY);
            }
        },
        {"TPCData"});

    auto* c0 {new TCanvas {"c0", "Bragg peak"}};
    hBragg->Draw("hist");
    auto* c1 = new TCanvas("c1", "beam");
    h2d->Draw("colz");


    if(compare_with_simu)
    {

        auto* srim {new ActPhysics::SRIM()};
        srim->ReadTable("beamInGas", "../../Simulation/SRIM/17F_H2-iC4H10_95-5_760mbar.txt");
        ActPhysics::Particle p {"17F"};

        std::vector<double> energies;
        auto hSim = new TH1D {"hSim", "", 128, 0, 128};
        auto dsim = gated
                        .Define("Tini",
                                [&](ROOT::Math::XYZPointF& lastPoint)
                                {
                                    float x {lastPoint.X()};
                                    auto Tini {srim->EvalInitialEnergy("beamInGas", 0, x)};
                                    return Tini;
                                },
                                {"lastPoint"})
                        .Define("lastPX", "lastPoint.X()");

        dsim.Foreach(
            [&](double Tini, float x)
            {
                double currentEne = Tini;
                for(int i=0; i < int(x); i++)
                {
                    double dEdx = srim->EvalStoppingPower("beamInGas", currentEne);
                    currentEne = srim->Slow("beamInGas", currentEne, 1);
                    hSim->Fill(i, dEdx);
                }
            },
            {"Tini", "lastPX"});

        // auto* c2 = new TCanvas("c2", "simu");
        hSim->Scale(hBragg->GetMaximum()/hSim->GetMaximum());
        hSim->SetLineColor(2);
        c0->cd();
        hSim->DrawClone("hist same");
    }
}
