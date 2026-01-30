// inspect potential candidate evts for p,2p reaction
// since MergerDetector only works for 2-body vertices, we need to do may of the calculations done by the MergerDetector
// class (TL, angles, etc.)

#include "ActDataManager.h"
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActOptions.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTypes.h"

#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include <fstream>

void ScalePoint(ROOT::Math::XYZPointF& point, float xy, float z) // void ActRoot::MergerDetector::ScalePoint function
                                                                 // from $ACTROOT/Detectors/src/ActMergerDetector.cxx
{
    point += ROOT::Math::XYZVectorF(0.5f, 0.5f, 0.5f);
    point.SetX(point.X() * xy);
    point.SetY(point.Y() * xy);
    point.SetZ(point.Z() * z);
}

double GetTheta3D(const ROOT::Math::XYZVectorF& beam,
                  const ROOT::Math::XYZVectorF& other) // void ActRoot::MergerDetector::GetTheta3D function
                                                       // from $ACTROOT/Detectors/src/ActMergerDetector.cxx
{
    auto dot {beam.Unit().Dot(other.Unit())};
    return TMath::ACos(dot) * TMath::RadToDeg();
}

double GetPhi3D(const ROOT::Math::XYZVectorF& beam,
                const ROOT::Math::XYZVectorF& other) // void ActRoot::MergerDetector::GetPhi3D function
                                                     // from $ACTROOT/Detectors/src/ActMergerDetector.cxx
{
    auto trackUnitary {other.Unit()};
    return TMath::ATan2(trackUnitary.Y(), trackUnitary.Z()) * TMath::RadToDeg();
}

int inspect_candidates()
{
    // Read in candidate evts
    std::ifstream fIn("./Outputs/candidates.dat");
    std::vector<std::pair<double, double>> candidates;
    double run, entry;
    while(fIn >> run >> entry)
    {
        candidates.emplace_back(run, entry);
    }

    // Filter the data for those evts -- there has to be a better way to do this
    ActRoot::DataManager dataman {"../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EMerge)};
    chain->AddFriend(chain2.get());
    auto chain3 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain3.get());

    // ROOT::EnableImplicitMT();
    ROOT::RDataFrame d {*chain};

    auto df {d.Filter(
        [&](ActRoot::MergerData& mer)
        {
            for(const auto& [r, e] : candidates)
            {
                if(mer.fRun == r && mer.fEntry == e)
                    return true;
            }
            return false;
        },
        {"MergerData"})};


    // grab drift factor, we'll need it for conversions
    ActRoot::InputParser parser {"../configs/detector.conf"};
    auto bl {parser.GetBlock("Merger")};
    auto drift {bl->GetDouble("DriftFactor")};
    std::cout << drift << std::endl;


    // Express track length as energy with SRIM
    // auto* srim {new ActPhysics::SRIM()};
    // srim->ReadTable("alphasInGas", "../../Simulation/SRIM/4He_H2-iC4H10_95-5_755mbar.txt");

    // Calculate quantities the MergerDetector would for binary reactions
    df = df.Define("Qtot",
                   [](const ActRoot::TPCData& tpc)
                   {
                       std::vector<double> qtot;
                       qtot.reserve(tpc.fClusters.size());

                       for(const auto& cl : tpc.fClusters)
                       {
                           //    if(cl.GetIsBeamLike())
                           //        continue;
                           auto q {0.};
                           for(const auto& vx : cl.GetVoxels())
                               q += vx.GetCharge();

                           qtot.push_back(q);
                       }
                       return qtot;
                   },
                   {"TPCData"})
             .Define("EndPoint",
                     [&](const ActRoot::TPCData& tpc)
                     {
                         std::vector<ROOT::Math::XYZPointF> end_v;
                         end_v.reserve(tpc.fClusters.size());
                         for(auto& cl : tpc.fClusters)
                         {
                             //  if(cl.GetIsBeamLike())
                             //      continue;
                             auto cl_copy = cl; // make a copy because we will modify (sort) and cl is inherited as
                                                // const from const variable tpc
                             auto line = cl_copy.GetRefToLine();
                             auto dir = line.GetDirection();
                             cl_copy.SortAlongDir(dir);
                             auto end = cl_copy.GetRefToVoxels().back().GetPosition();
                             ScalePoint(end, 2., drift);
                             line.Scale(2, drift);
                             auto projEnd = line.ProjectionPointOnLine(end);
                             end_v.push_back(projEnd);
                         }
                         return end_v;
                     },
                     {"TPCData"})
             .Define("TL", // find TL -- similar to ActRoot::MergerDetector::TrackLengthFromLightIt
                     [&](const ActRoot::TPCData& tpc)
                     {
                         std::vector<double> tl;
                         tl.reserve(tpc.fClusters.size());
                         for(auto& cl : tpc.fClusters)
                         {
                             if(cl.GetIsBeamLike())
                                 tl.push_back(0.);
                             else
                             {
                                 auto cl_copy = cl;
                                 auto begin = tpc.fRPs[0];
                                 ScalePoint(begin, 2., drift);
                                 auto line = cl_copy.GetRefToLine();
                                 auto dir = line.GetDirection();
                                 cl_copy.SortAlongDir(dir);
                                 auto end = cl_copy.GetRefToVoxels().back().GetPosition();
                                 ScalePoint(end, 2., drift);
                                 line.Scale(2, drift);
                                 auto projBegin = line.ProjectionPointOnLine(begin);
                                 auto projEnd = line.ProjectionPointOnLine(end);
                                 tl.push_back((projBegin - projEnd).R());
                             }
                         }
                         return tl;
                     },
                     {"TPCData"})
             .Define("Angles", // find angles theta and phi -- similar to ActRoot::MergerDetector::ComputeAngles
                     [&](const ActRoot::TPCData& tpc)
                     {
                         std::vector<std::pair<double, double>> angles; // theta, phi
                         angles.reserve(tpc.fClusters.size());
                         // find beam direction
                         ROOT::Math::XYZVectorF beamDir {};
                         for(auto& cl : tpc.fClusters)
                         {
                             if(cl.GetIsBeamLike())
                                 beamDir = cl.GetLine().GetDirection().Unit();
                         }
                         for(auto& cl : tpc.fClusters)
                         {
                             if(cl.GetIsBeamLike())
                                 angles.push_back({0., 0.});
                             else
                             {
                                 auto theta {GetTheta3D(beamDir, cl.GetLine().GetDirection())};
                                 auto phi {GetPhi3D(beamDir, cl.GetLine().GetDirection())};
                                 angles.push_back({theta, phi});
                             }
                         }
                         return angles;
                     },
                     {"TPCData"})
             .Define("IsRecoil",
                     [&](const ActRoot::TPCData& tpc, const std::vector<double> q)
                     {
                         std::vector<bool> isRecoil(tpc.fClusters.size(), false);
                         ROOT::Math::XYZVectorF beamDir {};
                         for(auto& cl : tpc.fClusters)
                         {
                             if(cl.GetIsBeamLike())
                                 beamDir = cl.GetLine().GetDirection().Unit();
                         }


                         int iq = 0;
                         int idxMaxQ = -1;
                         int idxMinAng = -1;
                         double maxQ = -1.;
                         double minAng = 1e9;

                         for(size_t i = 0; i < tpc.fClusters.size(); ++i)
                         {
                             const auto& cl = tpc.fClusters[i];
                             if(cl.GetIsBeamLike())
                                 continue;

                             // highest Qtot
                             if(iq < (int)q.size() && q[iq] > maxQ)
                             {
                                 maxQ = q[iq];
                                 idxMaxQ = i;
                             }

                             // smallest angle to beam
                             auto dir = cl.GetLine().GetDirection().Unit();
                             double ang = std::acos(dir.Dot(beamDir));
                             if(ang < minAng)
                             {
                                 minAng = ang;
                                 idxMinAng = i;
                             }

                             ++iq;
                         }
                         if(idxMaxQ == idxMinAng && idxMaxQ >= 0)
                         {
                             isRecoil[idxMaxQ] = true;
                         }
                         else
                         {
                             std::cerr << "I'm confused: maxQ idx=" << idxMaxQ << " minAng idx=" << idxMinAng
                                       << std::endl;
                         }

                         return isRecoil;
                     },
                     {"TPCData", "Qtot"});


    /*
    void ActRoot::MergerDetector::ComputeAngles()
    {

        XYZVector beamDir {};
        if(fBeamPtr)
            beamDir = fBeamPtr->GetLine().GetDirection().Unit();
        else
            beamDir = {1, 0, 0};
        // Light
        // 1-> Thetas (under different beam directions)
        fMergerData->fThetaLight = GetTheta3D(beamDir, fLightPtr->GetLine().GetDirection());
        fMergerData->fThetaLegacy =
            fMergerData->fThetaLight; // this is a copy just in case we apply the Corrector det (only for E796)
        // Debug: angle computed assuming beam exactly along X axis
        fMergerData->fThetaDebug = GetTheta3D({1, 0, 0}, fLightPtr->GetLine().GetDirection());
        // 2-> Phi
        fMergerData->fPhiLight = GetPhi3D(beamDir, fLightPtr->GetLine().GetDirection());

        // Beam
        fMergerData->fThetaBeam = GetTheta3D({1, 0, 0}, beamDir);
        fMergerData->fThetaBeamZ = TMath::ATan(beamDir.Z() / beamDir.X()) * TMath::RadToDeg();
        fMergerData->fPhiBeamY = TMath::ATan(beamDir.Y() / beamDir.X()) * TMath::RadToDeg();

        // Heavy
        if(fHeavyPtr)
        {
            fMergerData->fThetaHeavy = GetTheta3D(beamDir, fHeavyPtr->GetLine().GetDirection());
            fMergerData->fPhiHeavy = GetPhi3D(beamDir, fHeavyPtr->GetLine().GetDirection());
        }
    }
    */

    // auto dfEne = dfconstZ.Define("Ene",
    //                              [&](double Lxy)
    //                              {
    //                                  double eneMeV {srim->EvalInitialEnergy("alphasInGas", 0, Lxy)};
    df.Foreach(
        [&](ActRoot::TPCData& tpc, ActRoot::MergerData& mer, std::vector<double> qtot, std::vector<double> TL,
            std::vector<std::pair<double, double>> angles, std::vector<bool> isRecoil)
        {
            std::cout << "Event: " << mer.fRun << " " << mer.fEntry << std::endl;
            std::cout << "N of clusters: " << tpc.fClusters.size() << std::endl;
            std::cout << "N of RPs: " << tpc.fRPs.size() << std::endl;
            auto begin {tpc.fRPs[0]};
            ScalePoint(begin, 2., drift);
            std::cout << "RP: " << tpc.fRPs[0] << " or scaled " << begin << std::endl;
            std::cout << "Ejectiles: " << std::endl;


            for(int i = 0; i < tpc.fClusters.size(); i++)
            {
                auto cl {tpc.fClusters[i]};
                // if(!cl.GetIsBeamLike())
                // {

                // // find TL -- similar to ActRoot::MergerDetector::TrackLengthFromLightIt
                // auto begin {tpc.fRPs[0]};
                // ScalePoint(begin, 2., drift);
                auto line {cl.GetRefToLine()};
                auto dir {line.GetDirection()};
                cl.SortAlongDir(dir);
                auto end {cl.GetRefToVoxels().back().GetPosition()};
                ScalePoint(end, 2., drift);
                line.Scale(2, drift);
                // auto projBegin {line.ProjectionPointOnLine(begin)};
                auto projEnd {line.ProjectionPointOnLine(end)};
                // auto TL {(projBegin - projEnd).R()};


                // std::cout << "start: "<<projBegin<<", end: "<<projEnd<<", TL: "<<TL << std::endl;
                if(cl.GetIsBeamLike())
                    std::cout << std::setprecision(4) << i << ": Beam, \tend: " << projEnd << ",\tQtot: " << qtot[i]
                              << std::endl;
                else if(isRecoil[i])
                    std::cout << std::setprecision(4) << i << ": Recoil, \tend: " << projEnd << ",\tQtot: " << qtot[i]
                              << ", \tTL: " << TL[i] << ", \t(theta, phi): (" << angles[i].first << " ,"
                              << angles[i].second << ")" << std::endl;
                else
                    std::cout << std::setprecision(4) << i << ": Light, \tend: " << projEnd << ",\tQtot: " << qtot[i]
                              << ", \tTL: " << TL[i] << ", \t(theta, phi): (" << angles[i].first << " ,"
                              << angles[i].second << ")" << std::endl;
                // }
            }
            std::cout << "-----------------------------------------" << std::endl;
        },
        {"TPCData", "MergerData", "Qtot", "TL", "Angles", "IsRecoil"});

    return 1;
}