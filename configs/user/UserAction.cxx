#include "UserAction.h"

#include "ActColors.h"
#include "ActInputParser.h"
#include "ActTPCData.h"

#include "TMath.h"

#include <functional>
#include <memory>
#include <set>
#include <vector>

void ActAlgorithm::UserAction::UserAction::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    if(block->CheckTokenExists("MaxAngle"))
        fMaxAngle = block->GetDouble("MaxAngle");
    if(block->CheckTokenExists("MinLength"))
        fMinLength = block->GetDouble("MinLength");
    if(block->CheckTokenExists("CylinderR"))
        fCylinderR = block->GetDouble("CylinderR");
}

void ActAlgorithm::UserAction::Run()
{
    auto& clusters {fTPCData->fClusters};

    // Set storing clusters that are fixable
    std::set<int> references {};
    for(int i = 0; i < clusters.size(); i++)
    {
        auto& cluster {clusters[i]};
        auto [xmin, xmax] {cluster.GetXRange()};
        auto theta {TMath::ACos(cluster.GetLine().GetDirection().Unit().Dot(XYZVectorF {1, 0, 0}))};
        if((std::abs(theta) <= fMaxAngle * TMath::DegToRad()) && ((xmax - xmin) > fMinLength))
        {
            references.insert(i);
            if(fIsVerbose)
            {
                std::cout << BOLDGREEN << "-- UserAction: FixBeam -- " << '\n';
                std::cout << "  ref cluster " << i << '\n';
                std::cout << "  theta  : " << theta * TMath::RadToDeg() << '\n';
                std::cout << "  length : " << (xmax - xmin) << RESET << '\n';
            }
        }
    }
    // And now iterate over references and try to add all other clusters to it
    std::set<int, std::greater<int>> toDelete;
    for(const auto& ref : references)
    {
        auto& refcluster {clusters[ref]};

        for(int i = 0; i < clusters.size(); i++)
        {
            if(i == ref) // do not compare equal clusters
                continue;
            if(toDelete.find(i) != toDelete.end()) // cluster already merged
                continue;
            // Get iter cluster
            auto& itcluster {clusters[i]};
            // Algorithm: get fit of reference cluster
            auto& line {refcluster.GetRefToLine()};
            // And check if center of our cluster is in fCylinder distance
            auto point {itcluster.GetLine().GetPoint()};
            auto dist {line.DistanceLineToPoint(point)};
            if(dist < fCylinderR)
            {
                if(fIsVerbose)
                {
                    std::cout << BOLDGREEN << "-- UserAction: FixBeam -- " << '\n';
                    std::cout << "  it cluster " << i << " is being merged with ref " << ref << '\n';
                    std::cout << "  dist : " << dist << RESET << '\n';
                }
                // Merge!!
                auto& saveVoxels {refcluster.GetRefToVoxels()};
                auto& delVoxels {itcluster.GetRefToVoxels()};
                saveVoxels.insert(saveVoxels.end(), std::make_move_iterator(delVoxels.begin()),
                                  std::make_move_iterator(delVoxels.end()));
                // Refit and recompute ranges
                refcluster.ReFit();
                refcluster.ReFillSets();
                // Mark to delete afterwards!
                toDelete.insert(i);
            }
        }
    }
    // Delete!
    for(const auto& idx : toDelete) // toDelete is sorted in greater order
        clusters.erase(clusters.begin() + idx);
}

void ActAlgorithm::UserAction::Print() const
{
    std::cout << BOLDCYAN << "····· " << GetActionID() << " ·····" << '\n';
    if(!fIsEnabled)
    {
        std::cout << "······························" << RESET << '\n';
        return;
    }
    std::cout << "  MaxAngle       : " << fMaxAngle << '\n';
    std::cout << "  MinLength      : " << fMinLength << '\n';
    std::cout << "  CylinderRadius : " << fCylinderR << RESET << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::UserAction* CreateUserAction()
{
    return new ActAlgorithm::UserAction;
}
