#include "L1Filter.h"

#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActTPCData.h"
#include "ActMergerData.h"
#include "ActMergerDetector.h"

#include <memory>

//
// AT 2026
// There are a lot of L1 events that are triggered by noise instead of a real light particle.
// Increase the maxVoxel variable of the [CleanDeltas] function can remove useful events with very small tracks,
// so we need to treat noise more carefully.
// This function acts after [FindRP].
//
// There are two types of events we want to tackle:
// 1. small clusters that are far in time from the RP. Reaction products should be correlated with the RP.
// 2. small clusters that are inside the beam path (under dev -- I ended up moving this part to Pipe0_PreProcess, cause I need the light particle ID)
//
void ActAlgorithm::L1Filter::L1Filter::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if (!fIsEnabled)
        return;
    if (block->CheckTokenExists("TLmax"))
        fTLmax = block->GetDouble("TLmax");
    if (block->CheckTokenExists("dZmin"))
        fdZmin = block->GetDouble("dZmin");
    // if (block->CheckTokenExists("OutOfConeVoxels"))
    //     fOutOfConeVoxels = block->GetDouble("OutOfConeVoxels");
    // if (block->CheckTokenExists("BeamMaxSpread"))
    //     fBeamMaxSpread = block->GetDouble("BeamMaxSpread");
}

// ActAlgorithm::L1Filter::BeamCone
// ActAlgorithm::L1Filter::GetBeamCone(const ActRoot::Cluster &cl)
// {
//     BeamCone cone;

//     ActRoot::Cluster c = cl;
//     c.SortAlongDir();

//     int nBins = 12;

//     cone.x0 = c.GetVoxels().front().GetPosition().X();
//     cone.xMax = c.GetVoxels().back().GetPosition().X();
//     cone.dxStep = (cone.xMax - cone.x0) / nBins;
//     cone.line = c.GetRefToLine(); // cone center axis is the beam axis
//     cone.yzWidth.resize(nBins, {0.0, 0.0});

//     for (const auto &v : c.GetVoxels())
//     {
//         double x = v.GetPosition().X();
//         double y = v.GetPosition().Y();
//         double z = v.GetPosition().Z();

//         int i = (x - cone.x0) / cone.dxStep;
//         if (i < 0 || i >= nBins)
//             continue;

//         auto lineProjection = cone.line.ProjectionPointOnLine(v.GetPosition()); // project voxel on beamline
//         double yLine = lineProjection.Y();
//         double zLine = lineProjection.Z();

//         // width of cone will be the max distance of each voxel from the
//         cone.yzWidth[i].first = std::max(cone.yzWidth[i].first, std::abs(y - yLine));
//         cone.yzWidth[i].second = std::max(cone.yzWidth[i].second, std::abs(z - zLine));
//     }

//     return cone;
// }

// bool ActAlgorithm::L1Filter::BeamCone::Contains(ActRoot::Voxel &v) const
// {
//     auto vPos = v.GetPosition();
//     int i = (vPos.X() - x0) / dxStep;
//     if (i < 0 || i >= (int)yzWidth.size())
//         return false;

//     auto [dy, dz] = yzWidth[i];
//     auto lineProjection = line.ProjectionPointOnLine(vPos);

//     return std::abs(vPos.Y() - lineProjection.Y()) <= dy &&
//            std::abs(vPos.Z() - lineProjection.Z()) <= dz;
// }

// double ActAlgorithm::L1Filter::BeamCone::GetDistfromCone(ActRoot::Voxel &v) const
// {
//     auto vPos = v.GetPosition();
//     int i = (vPos.X() - x0) / dxStep;
//     auto [dy, dz] = yzWidth[i];

//     auto lineProjection = line.ProjectionPointOnLine(vPos);

//     double yDistfromLine = std::abs(vPos.Y() - lineProjection.Y());
//     double zDistfromLine = std::abs(vPos.Z() - lineProjection.Z());

//     return std::max(yDistfromLine - dy, zDistfromLine - dz);
// }

void ActAlgorithm::L1Filter::Run()
{
    if (!fIsEnabled)
        return;

    if (fIsVerbose)
        std::cerr << BOLDMAGENTA << "=========================== L1Filter ===========================" << RESET << '\n';
    if (!fTPCData || fTPCData->fClusters.empty())
    {
        if (fIsVerbose)
            std::cerr << BOLDMAGENTA << "fTPCData is null or there are no clusters" << RESET << '\n';
        return;
    }

    // only look at L1 events
    if (fTPCData->fTrigger == 8) // assumes data is being processed with artroot -r tpc --with-trigger
    {
        // copy/ref clusters
        auto &clusters = fTPCData->fClusters;
        // get RP and beamlike particle
        if (fTPCData->fRPs.size() > 0)
        {
            // // calculate size of BL cone in case there's a lot of spread
            // int beamId = -1;
            // for (auto cl : clusters)
            // {
            //     if (cl.GetIsBeamLike())
            //         beamId = cl.GetClusterID();
            // }
            // auto beamCone = GetBeamCone(clusters[beamId]);

            for (const auto &rp : fTPCData->fRPs)
            {
                for (auto it = fTPCData->fClusters.begin(); it != fTPCData->fClusters.end();)
                {
                    // skip beamLike
                    if (it->GetIsBeamLike())
                    {
                        ++it;
                        continue;
                    }

                    it->SortAlongDir();
                    auto begin{it->GetVoxels().front().GetPosition()};
                    auto end{it->GetRefToVoxels().back().GetPosition()};

                    // // find how many voxels are out of the beam cone, and what's the maximum distance from the beam cone
                    // int voxelsOutOfCone{0};
                    // double maxDistFromCone{0};
                    // for (auto v : it->GetVoxels())
                    // {
                    //     if (!beamCone.Contains(v))
                    //     {
                    //         voxelsOutOfCone++;
                    //         maxDistFromCone = std::max(beamCone.GetDistfromCone(v), maxDistFromCone);
                    //     }
                    // }
                    // bool isInBeamCone{voxelsOutOfCone <= fOutOfConeVoxels};
                    // // bool isNearBeamCone{maxDistFromCone <= fBeamMaxSpread};
                    // bool isNearBeamCone{ !it->GetIsRecoil() && (maxDistFromCone <= fBeamMaxSpread)};

                    // find min dist from RP
                    auto distZ = TMath::Min(TMath::Abs(begin.Z() - rp.Z()), TMath::Abs(end.Z() - rp.Z()));
                    // and prelim TL (tracks I'm trying to filter are noise like)
                    auto tl{(begin - end).R()};

                    bool isShort{tl <= fTLmax};
                    bool isFarInZ{distZ >= fdZmin};

                    // Delete if conditions are acomplished
                    // if ((isShort && isFarInZ) || (isInBeamCone || isNearBeamCone))
                    if (isShort && isFarInZ)
                    {
                        it = fTPCData->fClusters.erase(it);
                        if (fIsVerbose)
                        {
                            std::cout << BOLDMAGENTA << "Deleted Cluster " << it->GetClusterID() << "\n";
                            std::cout << "TL: " << tl << "\n";
                            std::cout << "distZ: " << distZ << "\n";
                            // std::cout << "voxelsOutOfCone: " << voxelsOutOfCone << "\n";
                            // std::cout<<"Was it recoil ?" << " " <<"\n";
                            // std::cout << "maxDistFromCone: " << maxDistFromCone << RESET << "\n";
                        }
                    }
                    else
                    {
                        if (fIsVerbose)
                        {
                            std::cout << BOLDMAGENTA << "I didn't delete Cluster " << it->GetClusterID() << "\n";
                            std::cout << "TL: " << tl << "\n";
                            std::cout << "distZ: " << distZ << "\n";
                            // std::cout << "voxelsOutOfCone: " << voxelsOutOfCone << "\n";
                            // std::cout << "maxDistFromCone: " << maxDistFromCone << RESET << "\n";
                        }
                        ++it;
                    }
                }
            }
        }
    }

    if (fIsVerbose)
        std::cerr << BOLDMAGENTA << "================================================================" << RESET << '\n';
}

void ActAlgorithm::L1Filter::Print() const
{
    std::cout << BOLDMAGENTA << "····· " << GetActionID() << " ·····" << '\n';
    std::cout << " IsEnabled ? " << fIsEnabled << '\n';
    std::cout << " TLmax: " << fTLmax << '\n';
    std::cout << " dZmin: " << fdZmin << '\n';
    // std::cout << " OutOfConeVoxels: " << fOutOfConeVoxels << '\n';
    // std::cout << " BeamMaxSpread: " << fBeamMaxSpread << '\n';
    std::cout << "······························" << RESET << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::L1Filter *CreateUserAction()
{
    return new ActAlgorithm::L1Filter;
}
