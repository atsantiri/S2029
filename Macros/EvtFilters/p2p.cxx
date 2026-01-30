#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>

void p2p()
{
    // the three body vertices exist only in the filter data. So we start from there.
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    // The run numbers and event entry exist in the merger data so we want that too.
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    // We make the two chains friends so we can access the run and entry numbers from the filter data.
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());

    //// keep multipthreading commented out when writing in txt file
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};

    // Gate on candidates for (p,2p)
    auto df_filtered {
        df.Filter(
            [](ActRoot::TPCData& tpc, ActRoot::ModularData& m)
            {
                return (tpc.fClusters.size() == 4) && (tpc.fRPs.size() == 1) && (m.Get("GATCONF") == 8);
            }, // we want the cluster of 4 (1beam + 3products) with one reaction point && from ModularData L1Ok
               // trigger.
            {"TPCData", "ModularData"})
        .Filter(
            [](ActRoot::TPCData& tpc)
            {
                auto& rp {tpc.fRPs[0]};
                double min {40}; // modify based on resonance location
                double max {80};
                return (min <= rp.X()) && (rp.X() <= max);
        },
        {"TPCData"})
    };

    // Write to file: run and entry numbers that passed the cuts for 3-body vertices within min and max RPx
    std::ofstream streamer {"../Outputs/p2p_list.txt"};
    df_filtered.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    streamer.close();

    // ActPhysics::Particle pb {"17F"};
    // auto mbeam {pb.GetMass()};
    // ActPhysics::Particle pt {"p"};
    // auto mtarget {pt.GetMass()};
    // ActPhysics::Particle pl {"p"};

    // double EBeam {3.84}; // MeV/u from SRIM interpolation of exp. 17F range

    // // srim files
    // auto* srim {new ActPhysics::SRIM};
    // srim->ReadTable("beamInGas", "../../Simulation/SRIM/17F_H2-iC4H10_95-5_760mbar.txt");


    // auto def {df_filtered
    //               .Define("Line",
    //                       [](ActRoot::TPCData& tpc)
    //                       {
    //                           auto cl {tpc.fClusters[0]};
    //                           auto l {cl.GetRefToLine()};
    //                           return l;
    //                       },
    //                       {"TPCData"})
    //               .Define("RPx", [](ActRoot::TPCData& tpc) { return tpc.fRPs[0].X(); }, {"TPCData"})
    //               .Define("Ebeam",
    //                       [&](float rpx)
    //                       {
    //                           double Ebeam = srim->Slow("beamInGas", EBeam * pb.GetAMU(), rpx*2); // rpx comes from TPCData that hasn't done the unit conversion yet. RPx here is in pads not mm.
    //                           return Ebeam;
    //                       },
    //                       {"RPx"})
    //               .Define("ECM",
    //                       [&](double ebeam)
    //                       {
    //                           double ecm = pt.GetMass() / (pb.GetMass() + pt.GetMass()) * ebeam;
    //                           return ecm;
    //                       },
    //                       {"Ebeam"})
    //               .Define("ECN", [&](double ecm) { return ecm + 3.923; }, {"ECM"})
    //               };


    // // Plot histogram of RPx for event with 3 vertices. Is there some structure here we could use??
    // auto hRPx {def.Histo1D({"hRPx", "RP x histo;RP.X [pads];Counts", 256, 0, 256}, "RPx")};
    // auto hEcm {def.Histo1D({"hEcm", "Center of mass energy; ECM [MeV];Counts", 100, 0, 5}, "ECM")};
    // auto hEcn {def.Histo1D({"hEcn", "Compound nucleus Ex; Ex [MeV];Counts", 100, 3, 8}, "ECN")};

    // auto df_res{def.Filter([&](double ex){return (ex>=6. && ex<=6.3);},{"ECN"})};

    // // std::ofstream streamer {"../Outputs/p2p_resonance.txt"};
    // // df_res.Foreach([&](ActRoot::MergerData& mer) { mer.Stream(streamer); }, {"MergerData"});
    // // streamer.close();

    // // Draw
    // auto* c0 {new TCanvas {"c0", "(p,2p) canvas"}};
    // c0->DivideSquare(2);
    // c0->cd(1);
    // hRPx->DrawClone();
    // c0->cd(2);
    // hEcn->DrawClone();
    // c0->SaveAs("../Outputs/p2p_Ex_20Ne.png");
}
