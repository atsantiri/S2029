#include "L1Filter.h"

#include "ActCluster.h"
#include "ActColors.h"
#include "ActInputParser.h"
#include "ActMergerData.h"
#include "ActMergerDetector.h"
#include "ActTPCData.h"

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
// 2. small clusters that are inside the beam path (under dev -- I ended up moving this part to Pipe0_PreProcess, cause
// I need the light particle ID)
//
void ActAlgorithm::L1Filter::L1Filter::ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block)
{
    fIsEnabled = block->GetBool("IsEnabled");
    if(!fIsEnabled)
        return;
    if(block->CheckTokenExists("TLmax"))
        fTLmax = block->GetDouble("TLmax");
    if(block->CheckTokenExists("dZmin"))
        fdZmin = block->GetDouble("dZmin");
}

void ActAlgorithm::L1Filter::Run()
{
    if(!fIsEnabled)
        return;

    if(fIsVerbose)
        std::cerr << BOLDMAGENTA << "=========================== L1Filter ===========================" << RESET << '\n';
    if(!fTPCData || fTPCData->fClusters.empty())
    {
        if(fIsVerbose)
            std::cerr << BOLDMAGENTA << "fTPCData is null or there are no clusters" << RESET << '\n';
        return;
    }

    // only look at L1 events
    if(fTPCData->fTrigger == 8) // assumes data is being processed with artroot -r tpc --with-trigger
    {
        // copy/ref clusters
        auto& clusters = fTPCData->fClusters;
        // get RP and beamlike particle
        if(fTPCData->fRPs.size() > 0)
        {

            for(const auto& rp : fTPCData->fRPs)
            {
                for(auto it = fTPCData->fClusters.begin(); it != fTPCData->fClusters.end();)
                {
                    // skip beamLike
                    if(it->GetIsBeamLike())
                    {
                        ++it;
                        continue;
                    }

                    it->SortAlongDir();
                    auto begin {it->GetVoxels().front().GetPosition()};
                    auto end {it->GetRefToVoxels().back().GetPosition()};

                    // find min dist from RP
                    auto distZ = TMath::Min(TMath::Abs(begin.Z() - rp.Z()), TMath::Abs(end.Z() - rp.Z()));
                    // and prelim TL (tracks I'm trying to filter are noise like)
                    auto tl {(begin - end).R()};

                    bool isShort {tl <= fTLmax};
                    bool isFarInZ {distZ >= fdZmin};

                    // Delete if conditions are acomplished
                    if(isShort && isFarInZ)
                    {
                        it = fTPCData->fClusters.erase(it);
                        if(fIsVerbose)
                        {
                            std::cout << BOLDMAGENTA << "Deleted Cluster " << it->GetClusterID() << "\n";
                            std::cout << "TL: " << tl << "\n";
                            std::cout << "distZ: " << distZ << "\n";
                        }
                    }
                    else
                    {
                        if(fIsVerbose)
                        {
                            std::cout << BOLDMAGENTA << "I didn't delete Cluster " << it->GetClusterID() << "\n";
                            std::cout << "TL: " << tl << "\n";
                            std::cout << "distZ: " << distZ << "\n";
                        }
                        ++it;
                    }
                }
            }
        }
    }
    else
    {
        if(fIsVerbose)
            std::cout << BOLDMAGENTA << "Not L1 event" << RESET << '\n';
    }

    if(fIsVerbose)
        std::cerr << BOLDMAGENTA << "================================================================" << RESET << '\n';
}

void ActAlgorithm::L1Filter::Print() const
{
    std::cout << BOLDMAGENTA << "····· " << GetActionID() << " ·····" << '\n';
    std::cout << " IsEnabled ? " << fIsEnabled << '\n';
    std::cout << " TLmax: " << fTLmax << '\n';
    std::cout << " dZmin: " << fdZmin << '\n';
    std::cout << "······························" << RESET << '\n';
}

// Create symbol to load class from .so
extern "C" ActAlgorithm::L1Filter* CreateUserAction()
{
    return new ActAlgorithm::L1Filter;
}
