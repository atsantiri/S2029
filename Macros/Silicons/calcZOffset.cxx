#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilSpecs.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TF1.h"

#include <fstream>

void calcZOffset()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EMerge};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get());

    ActPhysics::SilSpecs specs;
    specs.ReadFile("../../configs/silspecs.conf");

    ROOT::RDataFrame d {*chain};

    auto df {d.Filter([&](ActRoot::TPCData& tpc, ActRoot::ModularData& mod)
                      { return (mod.Get("GATCONF") == 64) && (tpc.fClusters.size() == 1); },
                      {"TPCData", "ModularData"})};

    std::vector<double> zpos, ypos;
    auto* h2d {new TH2D {"hPad", "Pad plane;X [pad];Z [pad]", 128, 0, 128, 128, 0, 128}};
    auto* hZ {new TH1D {"hZ", "Z pos;Counts;Z [pad]", 100, 60, 128}};
    auto* hY {new TH1D {"hY", "Y pos;Counts;Y [pad]", 128, 0, 128}};
    df.Foreach(
        [&](ActRoot::TPCData& tpc)
        {
            auto& cluster {tpc.fClusters.front()};
            for(const auto& v : cluster.GetVoxels())
            {
                auto& pos {v.GetPosition()};
                if(pos.X() <= 10)
                {
                    zpos.push_back(pos.Z());
                    ypos.push_back(pos.Y());
                    hZ->Fill(pos.Z());
                    hY->Fill(pos.Y());
                    // h2d->Fill(pos.X(), pos.Z());}
                }
            }
        },
        {"TPCData"});

    double mean_zpos {TMath::Mean(zpos.size(), zpos.data())};
    std::cout << "average zpos: " << mean_zpos << std::endl;
    double mean_ypos {TMath::Mean(ypos.size(), ypos.data())};
    std::cout << "average ypos: " << mean_ypos << std::endl;

    // get fDriftFactor to convert pos.Z to mm
    ActRoot::InputParser parser {"../../configs/detector.conf"};
    auto block {parser.GetBlock("Merger")};
    auto fDriftFactor {block->GetDouble("DriftFactor")};
    auto mean_zpos_mm {mean_zpos * fDriftFactor};

    auto* c0 {new TCanvas {"c0", ""}};
    // h2d->Draw("colz");
    c0->DivideSquare(2);
    c0->cd(1);
    hZ->Draw();
    TF1* fz = new TF1("fz", "gaus(0)", 60, 120);
    fz->SetParameters(1.5E6, 84, 0.06);
    hZ->Fit(fz, "R");
    auto meanZ {fz->GetParameter(1)};
    std::cout << "mean Z pos: " << std::fixed << std::setprecision(3) << meanZ << " +- " << fz->GetParameter(2)
              << std::endl;
    std::cout<< "mean z pos in mm "<<mean_zpos_mm<<std::endl;
    c0->cd(2);
    hY->Draw();
    TF1* fy = new TF1("fy", "gaus(0)", 0, 128);
    // fy->SetParameters(1.5E6, 84, 0.06);
    hY->Fit(fy, "R");
    auto meanY {fy->GetParameter(1)};
    std::cout << "mean Y pos: " << std::fixed << std::setprecision(3) << meanY << " +- " << fy->GetParameter(2)
              << std::endl;

    for(const auto& layer : {"f0", "l0", "r0"})
    {
        auto z {specs.GetLayer(layer).GetPoint().Z()};
        std::cout << layer << ": " << z << ", zOffset = " << mean_zpos_mm - z << std::endl;
    }
}