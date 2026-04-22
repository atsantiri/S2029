#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <stdexcept>

void calcMaxZdriftL1()
{
    std::string configsdir{"./../configs/"};

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman{configsdir + "data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{datman.GetJoinedData()};
    auto chain2{datman.GetChain(ActRoot::ModeType::EMerge)};
    auto chain3{datman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4{datman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    chain->AddFriend(chain3.get(), "TPCData");
    chain->AddFriend(chain4.get(), "GETTree");

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d{*chain};

    // Get GATCONF values
    auto df{d.Define("GATCONF", [](ActRoot::ModularData &mod)
                     { return static_cast<int>(mod.fLeaves["GATCONF"]); },
                     {"ModularData"})};

    // Get drift factor
    ActRoot::InputParser parserDet{configsdir + "detector.conf"};
    auto bl{parserDet.GetBlock("Merger")};
    auto drift{bl->GetDouble("DriftFactor")};

    // For L1 events: find max z distance of light particle from rp
    auto def = df.Define("zDrift", [&](ActRoot::MergerData &m, ActRoot::TPCData &tpc, int &gatconf)
                         {
                            
                            if (gatconf!=8)
                                return -300.0f;

                                // get light particle
                            auto idx = m.fLightIdx;
                            if (idx < 0)
                                return -300.0f;

                            if (static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                                return -300.0f;

                            auto &cluster = tpc.fClusters.at(idx);
                            auto &voxels = cluster.GetRefToVoxels();

                            // get RP
                            if (tpc.fRPs.empty())
                                return -300.0f;

                            auto rpVox = tpc.fRPs.front();

                            float maxDist = -1.0;
                            float zExtreme = 0.0;

                            for (auto &v : voxels)
                            {
                                auto pos = v.GetPosition();

                                float dx = pos.X() - rpVox.X();
                                float dy = pos.Y() - rpVox.Y();
                                float dz = pos.Z() - rpVox.Z();

                                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                                if (dist > maxDist)
                                {
                                    maxDist = dist;
                                    zExtreme = pos.Z();
                                }
                            }
                            float deltaZ = (zExtreme - rpVox.Z()) * drift;
                            float zDrift = deltaZ;
                            return zDrift; }, {"MergerData", "TPCData", "GATCONF"})
                   .Define("RunNumber",
                           [](ActRoot::MergerData &m) -> int
                           { return m.fRun; },
                           {"MergerData"});

    // 2D histo: x = run number, y = zDrift — bin width 1 in run so each bin = one run
    auto hZdriftVsRun = def.Histo2D(
        {"hZdriftVsRun", "Z drift vs Run;Run number;Z drift [mm]", 25, 45, 70, 400, -250, 250},
        "RunNumber", "zDrift");

    auto hZdrift = def.Histo1D({"hZdrift", "Z drift distribution;Z drift [mm];Counts", 400, -250, 250}, "zDrift");
    auto c = new TCanvas("c", "c", 800, 600);

    hZdrift->DrawClone();
    auto *cZdriftRun = new TCanvas("cZdriftRun", "Z drift vs Run (2D)", 800, 600);
    hZdriftVsRun->DrawClone("COLZ");
}
