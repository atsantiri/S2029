#include "ActDataManager.h"
#include "ActModularData.h"
#include "ActCutsManager.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TROOT.h"

#include <atomic>
#include <stdexcept>

// PreProcessing for S2029
// First step is to find events that hit the pad plane, by calculating the dZ distribution of the voxel with highest Z from the RP and removing the edge peak
// Second step is to try and remove the events triggered on noise caused from the spread of the beam. This can't be done by increasing maxVoxeln cause I'll lose tiny alpha tracks

struct BeamCone
{
    std::vector<std::pair<double, double>> yzWidth;
    double x0, xMax, dxStep;
    ActRoot::Line line;
    bool Contains(ActRoot::Voxel &v) const;
    double GetDistfromCone(ActRoot::Voxel &v) const;
};

BeamCone GetBeamCone(const ActRoot::Cluster &cl)
{
    BeamCone cone;

    ActRoot::Cluster c = cl;
    c.SortAlongDir();

    int nBins = 12;

    cone.x0 = c.GetVoxels().front().GetPosition().X();
    cone.xMax = c.GetVoxels().back().GetPosition().X();
    cone.dxStep = (cone.xMax - cone.x0) / nBins;
    cone.line = c.GetRefToLine(); // cone center axis is the beam axis
    cone.yzWidth.resize(nBins, {0.0, 0.0});

    for (const auto &v : c.GetVoxels())
    {
        double x = v.GetPosition().X();
        double y = v.GetPosition().Y();
        double z = v.GetPosition().Z();

        int i = (x - cone.x0) / cone.dxStep;
        if (i < 0 || i >= nBins)
            continue;

        auto lineProjection = cone.line.ProjectionPointOnLine(v.GetPosition()); // project voxel on beamline
        double yLine = lineProjection.Y();
        double zLine = lineProjection.Z();

        // width of cone will be the max distance of each voxel from the
        cone.yzWidth[i].first = std::max(cone.yzWidth[i].first, std::abs(y - yLine));
        cone.yzWidth[i].second = std::max(cone.yzWidth[i].second, std::abs(z - zLine));
    }

    return cone;
}

bool BeamCone::Contains(ActRoot::Voxel &v) const
{
    auto vPos = v.GetPosition();
    int i = (vPos.X() - x0) / dxStep;
    if (i < 0 || i >= (int)yzWidth.size())
        return false;

    auto [dy, dz] = yzWidth[i];
    auto lineProjection = line.ProjectionPointOnLine(vPos);

    return std::abs(vPos.Y() - lineProjection.Y()) <= dy &&
           std::abs(vPos.Z() - lineProjection.Z()) <= dz;
}

double BeamCone::GetDistfromCone(ActRoot::Voxel &v) const
{
    auto vPos = v.GetPosition();
    int i = (vPos.X() - x0) / dxStep;
    auto [dy, dz] = yzWidth[i];

    auto lineProjection = line.ProjectionPointOnLine(vPos);

    double yDistfromLine = std::abs(vPos.Y() - lineProjection.Y());
    double zDistfromLine = std::abs(vPos.Z() - lineProjection.Z());

    return std::max(yDistfromLine - dy, zDistfromLine - dz);
}

void Pipe0_PreProcess(const std::string &beam = "17F")
{

    //===========================================================
    bool filterZ = true;
    //===========================================================

    std::string configsdir{"./../configs/"};

    ROOT::EnableImplicitMT();
    // Read data
    ActRoot::DataManager datman{configsdir + "data.conf", ActRoot::ModeType::EReadSilMod};
    auto chain{datman.GetJoinedData()};
    auto chain2{datman.GetChain(ActRoot::ModeType::EMerge)};
    auto chain3{datman.GetChain(ActRoot::ModeType::EFilter)};
    auto chain4{datman.GetChain(ActRoot::ModeType::EReadTPC)};
    chain->AddFriend(chain2.get());
    // chain->AddFriend(chain3.get(), "TPCData");
    // chain->AddFriend(chain4.get(), "GETTree");
    chain->AddFriend(chain3.get());
    chain->AddFriend(chain4.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d{*chain};

    // Get drift factor
    ActRoot::InputParser parserDet{configsdir + "detector.conf"};
    auto bl{parserDet.GetBlock("Merger")};
    auto drift{bl->GetDouble("DriftFactor")};

    // For L1 events: find mean Z of beamlike particle
    auto df = d.Define("GATCONF", [](ActRoot::TPCData &tpc)
                       { return static_cast<int>(tpc.fTrigger); }, {"TPCData"})
                  .Define("zDrift", [&](ActRoot::MergerData &m, ActRoot::TPCData &tpc)
                          {
                            auto idx = m.fLightIdx;
                            if (idx < 0)
                                return -300.0f;

                            if (static_cast<std::size_t>(idx) >= tpc.fClusters.size())
                                return -300.0f;

                            auto &cluster = tpc.fClusters.at(idx);
                            auto &voxels = cluster.GetRefToVoxels();

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
                            return zDrift; }, {"MergerData", "TPCData"})
                  .Define("RunNumber",
                          [](ActRoot::MergerData &m) -> int
                          { return m.fRun; },
                          {"MergerData"});
                //   .Define("BeamCone", [&](ActRoot::TPCData &tpc)
                //           {
                //             int beamId = -1;
                //             for (auto cl : tpc.fClusters)
                //             {
                //                 if (cl.GetIsBeamLike())
                //                     beamId = cl.GetClusterID();
                //             }
                //             auto beamCone = GetBeamCone(tpc.fClusters[beamId]);
                //             return beamCone; }, {"TPCData"});


    auto hZdrift = df.Histo1D({"hZdrift", "Z drift distribution;Z drift [mm];Counts", 300, -200, 200}, "zDrift");
    auto c = new TCanvas("c", "c", 800, 600);
    hZdrift->DrawClone();

    // Results from fitting the two peaks near the edges of the zDrift distribution
    float maxZ_mean = 1.44452e+02;
    float maxZ_sig = 2.13408e+00;
    float minZ_mean = -112.752;
    float minZ_sig = 2.29358;

    // Filter out events with very high or low Z drift.
    auto dFiltered = df.Filter([&](float zdrift, int gatconf)
                               {
                                    if (gatconf != 8)
                                        return true;
                                    return zdrift >= (minZ_mean + 2 * minZ_sig) && zdrift <= (maxZ_mean - 2 * maxZ_sig); },
                               {"zDrift", "GATCONF"});

    auto hZdriftFiltered = dFiltered.Histo1D({"hZdriftFiltered", "Z drift distribution;Z drift [mm];Counts", 300, -200, 200}, "zDrift");
    hZdriftFiltered->SetLineColor(2);
    hZdriftFiltered->DrawClone("same");

    // Save filtered dataframe in output
    auto outName{TString::Format("./Outputs/tree_preProcessed_%s.root", beam.c_str())};
    if (filterZ)
    {
        std::cout << "Saving Z filtered PreProcessed_Tree in file : " << outName << '\n';
        dFiltered.Snapshot("PreProcessed_Tree", outName.Data());
    }
    else
    {
        std::cout << "Saving PreProcessed_Tree in file : " << outName << '\n';
        std::cout << "This data will not be Z filtered (bool flag set to false)" << std::endl;
        df.Snapshot("PreProcessed_Tree", outName.Data());
    }
}
