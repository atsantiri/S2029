#include "ActColors.h"
#include "ActCutsManager.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include "ROOT/TThreadedObject.hxx"
#include "Rtypes.h"
#include <random>

#include "TCanvas.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include "Math/Point3Dfwd.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

#include "../PostAnalysis/HistConfig.h"

void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 455;

    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

double AngleWithNormal(const ROOT::Math::XYZVector& dir, const ROOT::Math::XYZVector& normal)
{
    auto dot {dir.Unit().Dot(normal.Unit())};
    return TMath::ACos(dot);
}

std::pair<ROOT::Math::XYZPoint, ROOT::Math::XYZPoint>
SampleVertex(double meanZ, double sigmaZ, double meanY, double sigmaY, TH3D* h, double lengthX)
{

    // X is always common for both manners
    double Xstart {0};
    double Xrp {gRandom->Uniform() * lengthX};
    // Y depends completely on the method of calculation
    double Ystart {-1};
    double Yrp {-1};
    // Z of beam at entrance
    double Zstart {gRandom->Gaus(meanZ, sigmaZ)};
    double Zrp {-1};
    // Ystart in this case is sampled from the histogram itself!
    double thetaXY {};
    double thetaXZ {};

    // Two options depending on the existance of a emittance histogram or not
    if(h)
    {
        // Ystart in this case is sampled from the histogram itself!
        double thetaXY {};
        double thetaXZ {};
        h->GetRandom3(Ystart, thetaXY, thetaXZ);
        // Mind that Y is not centred in the histogram value!
        // Rp values are computed as follows:
        Yrp = Ystart - Xrp * TMath::Tan(thetaXY * TMath::DegToRad());
        Zrp = Zstart - Xrp * TMath::Tan(thetaXZ * TMath::DegToRad());
    }
    else
    {
        // Starting point in Y is also random
        Ystart = gRandom->Gaus(meanY, sigmaY);
        // And now rp values
        // This way, emittance is fully random
        Yrp = gRandom->Gaus(meanY, sigmaY);
        Zrp = gRandom->Gaus(meanZ, sigmaZ);
    }

    // Mind that Y is not centred in the histogram value!
    // Rp values are computed as follows:
    Yrp = Ystart - Xrp * TMath::Tan(thetaXY * TMath::DegToRad());
    Zrp = Zstart - Xrp * TMath::Tan(thetaXZ * TMath::DegToRad());
    ROOT::Math::XYZPoint start {Xstart, Ystart, Zstart};
    ROOT::Math::XYZPoint vertex {Xrp, Yrp, Zrp};
    return {std::move(start), std::move(vertex)};
}

void ApplyNaN(double& val, double thresh = 0, const std::string& comment = "stopped")
{
    if(val <= thresh)
        val = std::nan(comment.c_str());
}

void makeGrid(std::string layer, std::unordered_map<std::string, ActPhysics::SilMatrix*> smAll)
{
    auto& cuts = smAll[layer]->GetGraphs();
    for(const auto& [id, cut] : cuts)
    {
        if(cut)
            cut->Draw("same");
    }
}

void Simulation_S2029(const std::string& beam, const std::string& target, const std::string& light,
                      const std::string& heavy, double T1, double Ex, int pressure, bool standalone)
{
    // set batch mode if not an independent function
    if(!standalone)
        gROOT->SetBatch(true);

    double angle_min {130.};
    double angle_max {140.};
    double qval {3.923}; // qvalue of 17F+p

    TRandom random;

    // Resolutions
    const double sigmaSil {0.05}; // AT: To change, get sigma from calibrations
    const double sigmaPercentBeam {0};
    const double sigmaAngleLight {0.95 / 2.355};
    // Parameters of beam in mm
    // Beam has to be manually placed in the simulation
    // Centered in Z and Y with a width of 4 mm
    // Center in Z
    // AT: note that in simue756 zVertexMean is coming from silicon matrices
    const double zVertexMean {128. + 18.}; // beam not centered in chamber, upwards
                                           //    double zVertexMean {silCentre + beamOffset}; // in Miguel's

    const double zVertexSigma {4}; // AT, to get from beam
    // Center in Y - miguel's simu doesn't have Yvertex?
    const double yVertexMean {128.};
    const double yVertexSigma {4};
    // const double zVertexMean {83.59};
    // const double zVertexSigma {3.79};

    // THRESHOLDS FOR SILICONS -- AT double check from calibration files
    const double thresholdSi0 {1.};
    const double thresholdSi1 {1.};

    // number of iterations
    const int iterations {static_cast<int>(1e6)};

    // ACTIVATE STRAGGLING OR NOT
    bool stragglingInGas {true};
    bool stragglingInSil {true};
    bool silResolution {true};
    bool thetaResolution {true};

    // TPC basic parameters
    ActRoot::TPCParameters fTPC {"Actar"};

    // Init particles
    ActPhysics::Particle p1 {beam};
    ActPhysics::Particle p2 {target};
    ActPhysics::Particle p3 {light};
    ActPhysics::Particle p4 {heavy};

    // Init kinematics generator
    ActSim::KinematicGenerator kingen {p1, p2, p3, p4, 0, 0};
    kingen.Print();
    auto* kin {kingen.GetBinaryKinematics()};

    // Get threshold
    auto T1Thresh {ActPhysics::Kinematics(p1, p2, p3, p4, -1, Ex).GetT1Thresh()};


    // Silicon specs
    ActPhysics::SilSpecs specs;
    specs.ReadFile("../configs/silspecs.conf");

    // Silicon EFFECTIVE matrix
    double silCentre {};
    std::vector<std::string> silLayers {"f0", "l0", "r0"};
    TString secondLayer {"f1"};
    std::unordered_map<std::string, ActPhysics::SilMatrix*> smAll;
    for(const auto& l : silLayers)
        smAll[l] = specs.GetLayer(l).GetSilMatrix()->Clone();

    /////////////////////////////////////////////////////////////////////////////
    // need to move silicons around, like Ivan's lines 304 - 334 of s2384/Simulation/do_simu.cxx
    /////////////////////////////////////////////////////////////////////////////
    specs.DrawGeo();

    // // CUTS ON SILICON ENERGY, depending on particle and layer
    // not needed here, since only protons make it to the Silicons. If I had multiple particles in my pid this would be
    // needed. I may need it for L1
    // std::map<std::string, std::pair<double, double>> eLoss0Cuts {}; for(const auto& layer : {"f0", "l0",
    // "r0"})
    // {
    //     ActRoot::CutsManager<std::string> cut;
    //     cut.ReadCut(
    //         light,
    //         TString::Format("../PostAnalysis/Cuts/pid_%s_%s_%s.root", light.c_str(), layer, beam.c_str()).Data());
    //     if(cut.GetCut(light))
    //     {
    //         eLoss0Cuts[layer] = cut.GetXRange(light);
    //         std::cout << BOLDGREEN << "-> ESil in " << layer << " : " << light << ": [" << eLoss0Cuts[layer].first
    //                   << ", " << eLoss0Cuts[layer].second << "] MeV" << RESET << '\n';
    //     }
    //     else
    //     {
    //         std::cout << BOLDRED << "Simulation_S2029(): could not read PID cut for " << light
    //                   << " -> using default eLoss0Cut" << RESET << '\n';
    //         eLoss0Cuts[layer] = {0, 1000};
    //     }
    // }

    // Histograms
    // To compute a fine-grain efficiency, we require at least a binning width of 0.25 degrees!
    auto hThetaCM {HistConfig::ThetaCM.GetHistogram()};
    auto hThetaCMAll {HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaCM all", "All").GetHistogram()};
    auto hThetaLabAll {HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaLab all", "LabAll").GetHistogram()};
    auto hThetaLab {HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaLab", "Lab").GetHistogram()};
    auto hPhiAll {HistConfig::ChangeTitle(HistConfig::PhiCM, "PhiAll", "All").GetHistogram()};
    auto hPhiLab {HistConfig::ChangeTitle(HistConfig::PhiCM, "PhiLab", "Lab").GetHistogram()};

    auto hDistF0 {HistConfig::ChangeTitle(HistConfig::TL, "Distance to F0").GetHistogram()};
    auto hKinVertex {HistConfig::ChangeTitle(HistConfig::Kin, "Kinematics at vertex").GetHistogram()};
    auto hKinSampled {HistConfig::ChangeTitle(HistConfig::Kin, "Sampled kinematics").GetHistogram()};
    auto hEexBefore {HistConfig::ChangeTitle(HistConfig::Ex, "Ex before resolutions", "Bef").GetHistogram()};
    auto hEexAfter {HistConfig::ChangeTitle(HistConfig::Ex, "Ex after resolutions", "After").GetHistogram()};

    auto hRP {HistConfig::RP.GetHistogram()};
    auto hRPz {std::make_unique<TH2D>("hRPz", "RPz;Y [mm];Z [mm]", 550, 0, 256, 550, 0, 256)};
    auto hECN {HistConfig::ECN.GetHistogram()};
    auto hECNFr {HistConfig::ECN.GetHistogram()};  // front Sil
    auto hECNLat {HistConfig::ECN.GetHistogram()}; // lateral sil
    auto* hRPxSimu {HistConfig::RP.GetHistogram()->ProjectionX("hRPxSimu")};
    auto* hRPxRange {new TH2D {"hRPxRange", "Light range vs dist to sil wall;Dist to wall [mm];Light range [mm]", 300,
                               0, 500, 300, 0, 500}};
    auto* hEBeam {new TH2D {"hEBeam", "Range angle", 300, 0, 180, 300, 0, 500}};
    auto hRPxEBeam {HistConfig::EBeamRPx.GetHistogram()};

    auto hEffAll {HistConfig::ChangeTitle(HistConfig::Eff2D, "geometric 2D efficiency").GetHistogram()};
    auto hEffAfter {HistConfig::ChangeTitle(HistConfig::Eff2D, "Efficiency ~ E_{CM}").GetHistogram()};
    auto hEffDirAll {HistConfig::Eff2D.GetHistogram()};
    auto hEffDirAfter {HistConfig::Eff2D.GetHistogram()};
    hEffDirAfter->SetTitle("Efficiency;T_{Beam}^{Dir} [MeV];#theta_{3, Lab}^{dir} [#circ]");
    auto hEffECNAll {HistConfig::ChangeTitle(HistConfig::EcnThetaCM, "Efficiency ~ E_{^{18}Ne} All").GetHistogram()};
    auto hEffECNAfter {HistConfig::ChangeTitle(HistConfig::EcnThetaCM, "Efficiency ~ E_{^{18}Ne}").GetHistogram()};

    // ECM resoltion
    auto hECMRes {HistConfig::ECMECM.GetHistogram()};
    // T1Dir resolution
    auto hTBeamDirRes {HistConfig::ECMECM.GetHistogram()};
    hTBeamDirRes->SetTitle("T_{Beam} direct res;T_{Beam} sampled [MeV];T_{Beam} rec [MeV]");
    // auto hKin {HistConfig::Kin.GetHistogram()};
    // auto hECM {HistConfig::ECM.GetHistogram()};
    // auto* hprojection {new TH1D {"hprojection", "Ex;E_{x} [MeV]; Counts", 80, 11, 13}};
    // auto* htrack {
    //     new TH2D {"htrack", "trackslength;#theta_{heavy lab} [#circ];track length [mm]", 90, 0, 100, 90, 0, 250}};
    auto* hTestEcm {new TH2D {"hTestEcm", "E_{CM};E_{CM} ;E_{CM}Recon", 100, 0, 5, 100, 0, 5}};
    auto hTheta3CM {HistConfig::ThetaCM.GetHistogram()};
    auto hTheta3CMside {HistConfig::ThetaCM.GetHistogram()};
    auto hTheta3CMfront {HistConfig::ThetaCM.GetHistogram()};
    auto hTheta3Lab {HistConfig::ThetaCM.GetHistogram()};
    hTheta3Lab->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
    auto hTheta3Labside {HistConfig::ThetaCM.GetHistogram()};
    hTheta3Labside->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
    auto hTheta3Labfront {HistConfig::ThetaCM.GetHistogram()};
    hTheta3Labfront->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");

    TF1* fit = new TF1("fit", "[0]*(x^[1])", 10, -0.5); // fitting function for pid

    std::map<std::string, ROOT::TThreadedObject<TH2D>> hSilPID, hSilSP;
    auto hSP {HistConfig::SP.GetHistogram()};
    auto* hpid {new TH2D {"hPID", "PID;E[MeV];Eloss [MeV]", 1000, 0, 15, 1000, 0, 0.01}};
    for(const auto& layer : {"f0", "f1", "l0", "r0"})
    {
        hSilPID.emplace(layer, *hpid);
        hSilPID[layer]->SetTitle(TString::Format("%s", layer));
        hSilSP.emplace(layer, *hSP);
        hSilSP[layer]->SetTitle(TString::Format("%s", layer));
    }

    // Load SRIM tables
    // The name of the file sets particle + medium
    auto* srim {new ActPhysics::SRIM()};
    srim->ReadTable("light", TString::Format("SRIM/%s_H2-iC4H10_95-5_%dmbar.txt", light.c_str(), pressure).Data());
    srim->ReadTable("beam", TString::Format("SRIM/%s_H2-iC4H10_95-5_%dmbar.txt", beam.c_str(), pressure).Data());
    srim->ReadTable("lightInSil", TString::Format("SRIM/%s_silicon.txt", light.c_str()).Data());

    // Random generator
    auto* rand {new TRandom3()};
    rand->SetSeed(); // random path in each execution of macro

    // Runner: contains utility funcstions to execute multiple actions
    ActSim::Runner runner(nullptr, nullptr, rand, sigmaSil);

    // Output from simulation!
    // We only store a few things in the TTree
    // 1-> Excitation energy
    // 2-> Theta in CM frame
    // 3-> Weight of the generator: for three-body reactions (phase spaces) the
    // other two variables need to be weighted by this value. For binary
    // reactions, weight = 1
    // 4-> Energy at vertex
    // 5-> Theta in Lab frame
    auto* outFile {new TFile(TString::Format("Outputs/Simu_17F_p_p_%dmbar.root", pressure), "recreate")};
    auto* outTree {new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};

    double theta3CM_tree {};
    outTree->Branch("theta3CM", &theta3CM_tree);
    double Eex_tree {};
    outTree->Branch("Eex", &Eex_tree);
    double weight_tree {};
    outTree->Branch("weight", &weight_tree);
    double EVertex_tree {};
    outTree->Branch("EVertex", &EVertex_tree);
    double theta3Lab_tree {};
    outTree->Branch("theta3Lab", &theta3Lab_tree);
    double rpx_tree {};
    outTree->Branch("RPx", &rpx_tree);
    double ecn_tree {};
    outTree->Branch("ECN", &ecn_tree);

    //---- SIMULATION STARTS HERE
    ROOT::EnableImplicitMT();

    // timer
    TStopwatch timer {};
    timer.Start();

    // print fancy info
    std::cout << BOLDMAGENTA << "Running for Ex = " << Ex << " MeV" << RESET << '\n';
    std::cout << BOLDMAGENTA << "Output File: " << outFile->GetName() << RESET << '\n';
    std::cout << BOLDGREEN;
    const int percentPrint {5};
    int step {iterations / (100 / percentPrint)};
    int nextPrint {step};
    int percent {};
    double maxBinContent = -1;
    for(long int reaction = 0; reaction < iterations; reaction++)
    {
        // Print progress
        if(reaction >= nextPrint)
        {
            percent = 100 * (reaction + 1) / iterations;
            int nchar {percent / percentPrint};
            std::cout << "\r" << std::string((int)(percent / percentPrint), '|') << percent << "%";
            std::cout.flush();
            nextPrint += step;
        }
        // 1-> Sample vertex
        auto [start, vertex] {SampleVertex(zVertexMean, zVertexSigma, yVertexMean, yVertexSigma, nullptr, fTPC.X())};

        // 2-> Beam energy according to its sigma
        auto TBeam {runner.RandomizeBeamEnergy(
            T1 * p1.GetAMU(),
            sigmaPercentBeam * T1 * p1.GetAMU())}; // T1 in Mev / u * mass of beam in u = total kinetic energy

        // Slow down it according to vertex position
        auto distToVertex {(vertex - start).R()};
        // If distance is too small, protect from spline faults
        srim->SetUseSpline(false);
        if(stragglingInGas)
            TBeam = srim->SlowWithStraggling("beam", TBeam, distToVertex);
        else
            TBeam = srim->Slow("beam", TBeam, distToVertex);
        srim->SetUseSpline(true);

        // runner energy functions return std::nan when the particle is stopped in
        // the gas! if nan (aka stopped in gas, continue) if not stopped but beam
        // energy below kinematic threshold, continue
        if(std::isnan(TBeam) || TBeam < T1Thresh)
            continue;
        hRPxEBeam->Fill(vertex.X(), TBeam);

        ////////////////////////////////////////////////////////////////////////////////
        // 3-> Run kinematics!
        kingen.SetBeamAndExEnergies(TBeam, Ex);
        double theta3Lab {};
        double phi3Lab {};
        double T3Lab {};
        double theta4Lab {};
        double phi4Lab {};
        double weight {1};
        // Uniform phi always and it is the same for CM and Lab
        auto phiCM {gRandom->Uniform(0, TMath::TwoPi())};
        // thetaCM following xs or not
        double thetaCM = TMath::ACos(gRandom->Uniform(-1, 1));
        kin->ComputeRecoilKinematics(thetaCM, phiCM);

        // Light
        theta3Lab = kin->GetTheta3Lab();
        phi3Lab = phiCM;
        T3Lab = kin->GetT3Lab();
        // Heavy
        theta4Lab = kin->GetTheta4Lab();
        phi4Lab = phiCM;

        // Save without resolution
        // Simulated thetaCM, to be used for efficiency computation by CM interval and with our set reference direction
        double thetaCMEff {kin->ReconstructTheta3CMFromLab(T3Lab, theta3Lab)};
        double theta3LabEff {theta3Lab}; // before implementing resolution in angle
        // auto EexEff {kin->ReconstructExcitationEnergy(T3Lab, theta3Lab)};
        auto ECM {(p2.GetAMU() / (p2.GetAMU() + p1.GetAMU())) * TBeam};
        auto ECN {ECM + qval};
        // Beam energy in direct reaction
        double TBeamDir {kin->ComputeEquivalentOtherT1(TBeam)};
        // Theta3 in direct reaction
        double theta3LabDir {kin->ComputeOtherInLab(thetaCMEff).first * TMath::RadToDeg()};

        // Fill Histograms
        hThetaCMAll->Fill(thetaCMEff * TMath::RadToDeg());
        hThetaLabAll->Fill(theta3LabEff * TMath::RadToDeg());
        hPhiAll->Fill(phi3Lab * TMath::RadToDeg());
        hEffDirAll->Fill(theta3LabDir, TBeamDir);
        hEffAll->Fill(thetaCMEff * TMath::RadToDeg(), ECM);
        hRP->Fill(vertex.X(), vertex.Y());
        hEffECNAll->Fill(thetaCMEff * TMath::RadToDeg(),ECN);

        ////////////////////////////////////////////////////////////////////////////////
        // 4-> Include thetaLab resolution to compute thetaCM and Ex afterwards
        if(thetaResolution) // resolution in
            theta3Lab = gRandom->Gaus(theta3Lab, sigmaAngleLight * TMath::DegToRad());

        // Eval range of light particle
        auto lightRange {srim->EvalDirect("light", T3Lab)};
        // std::cout<<"lightRange = "<<lightRange<<std::endl;
        // auto heavyRange{srim->EvalDirect("heavy", T4Lab)};

        ////////////////////////////////////////////////////////////////////////////////
        // 5-> Propagate track from vertex to silicon wall using SilSpecs class
        // And using the angle with the uncertainty already in
        ROOT::Math::XYZVector dirBeamFrame {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                                            TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};

        // Declare beam direction
        auto beamDir {(vertex - start).Unit()};
        // Rotate to world = geometry frame
        auto dirWorldFrame {runner.RotateToWorldFrame(dirBeamFrame, beamDir)};

        // // Threshold L1, particles that stop in actar. Check before doing the continues
        // ROOT::Math::XYZPoint finalPointGas {vertex + lightRange * dirWorldFrame.Unit()};
        // if(0 <= finalPointGas.X() && finalPointGas.X() <= fTPC.X() && 0 <= finalPointGas.Y() &&
        //    finalPointGas.Y() <= fTPC.Y() && 0 <= finalPointGas.Z() && finalPointGas.Z() <= fTPC.Z())
        // {
        //     hSPCut->Fill(finalPointGas.Y(), finalPointGas.Z());
        //     // do something here????
        // }

        // Light Particle in Silicon
        int silIndex0 = -1;
        ROOT::Math::XYZPoint silPoint0; // in mm
        std::string layer0;
        std::shared_ptr<ActPhysics::SilMatrix> sm {};

        for(auto layer : silLayers)
        {
            auto [index, sp] = specs.FindSPInLayer(layer, vertex, dirWorldFrame);
            if(index != -1)
            {
                silIndex0 = index;
                silPoint0 = sp;
                layer0 = layer;
                sm = specs.GetLayer(layer).GetSilMatrix();
                break;
            }
        }

        // skip tracks that doesn't reach silicons or are in silicon index cut
        if(silIndex0 == -1)
            continue;

        // Apply SilMatrix cut
        if(TString(layer0).Contains("f"))
        {
            if(!sm->IsInside(silIndex0, silPoint0.Y(), silPoint0.Z()))
                continue;
        }
        else
        {
            if(!sm->IsInside(silIndex0, silPoint0.X(), silPoint0.Z()))
                continue;
        }

        // Define SP distance
        auto distance0 {(vertex - silPoint0).R()};
        // std::cout << "lightRange: " << lightRange << " distance0: " << distance0 << std::endl;
        double T3EnteringSil {-1};
        if(stragglingInGas)
            T3EnteringSil = srim->SlowWithStraggling("light", T3Lab, distance0);
        else
            T3EnteringSil = srim->Slow("light", T3Lab, distance0);
        // std::cout<<" SilPoint "<<silPoint0InMM<<", energy entering Sil: "<<T3EnteringSil<<std::endl;
        ApplyNaN(T3EnteringSil);
        // nan if stopped in gas
        if(!std::isfinite(T3EnteringSil))
            continue;

        // SILICON0
        ROOT::Math::XYZVector normal {specs.GetLayer(layer0).GetNormal()};
        auto angleNormal0 {AngleWithNormal(dirWorldFrame, normal)};
        double T3AfterSil0 {-1};
        if(stragglingInSil)
            T3AfterSil0 = srim->SlowWithStraggling("lightInSil", T3EnteringSil,
                                                   specs.GetLayer(layer0).GetUnit().GetThickness(), angleNormal0);
        else
            T3AfterSil0 =
                srim->Slow("lightInSil", T3EnteringSil, specs.GetLayer(layer0).GetUnit().GetThickness(), angleNormal0);

        auto eLoss0 {T3EnteringSil - T3AfterSil0};

        // Apply resolution
        if(T3AfterSil0 != 0)
        {
            if(silResolution)
                eLoss0 = gRandom->Gaus(eLoss0, sigmaSil * TMath::Sqrt(eLoss0 / 5.5));
            T3AfterSil0 = T3EnteringSil - eLoss0;
        }
        // ApplyNaN(eLoss0, thresholdSi0, "thresh");
        ApplyNaN(eLoss0, specs.GetLayer(layer0).GetThresholds().at(silIndex0));

        // nan if bellow threshold
        if(!std::isfinite(eLoss0))
            continue;

        ////////////////////////////////////////////////////////////////////////////////
        // 6-> Same but to silicon layer 1 if exists
        double T3AfterInterGas {};
        double distance1 {};
        int silIndex1 {};
        ROOT::Math::XYZPoint silPoint1 {};
        double eLoss1 {};
        double T3AfterSil1 {};
        bool isPunch {};
        if(T3AfterSil0 > 0 && (layer0 == "f0"))
        {
            // first, propagate in gas
            auto [silIndex1, silPoint1] {specs.FindSPInLayer(secondLayer.Data(), vertex, dirWorldFrame)};
            if(silIndex1 == -1)
                continue;

            distance1 = (silPoint1 - silPoint0).R();
            if(stragglingInSil)
                T3AfterInterGas = srim->SlowWithStraggling("light", T3AfterSil0, distance1);
            else
                T3AfterInterGas = srim->Slow("light", T3AfterSil0, distance1);
            ApplyNaN(T3AfterInterGas);
            if(!std::isfinite(T3AfterInterGas))
                continue;

            // now, second silicon layer if we have energy left
            if(T3AfterInterGas > 0)
            {
                // For S2008/S2029 angleNormal0 = angleNormal1 but this is not general
                auto angleNormal1 {angleNormal0};
                if(stragglingInSil)
                    T3AfterSil1 = srim->SlowWithStraggling("lightInSil", T3AfterInterGas,
                                                           specs.GetLayer(secondLayer.Data()).GetUnit().GetThickness(),
                                                           angleNormal1);
                else
                    T3AfterSil1 = srim->Slow("lightInSil", T3AfterInterGas,
                                             specs.GetLayer(secondLayer.Data()).GetUnit().GetThickness(), angleNormal1);

                auto eLoss1 {T3AfterInterGas - T3AfterSil1};
                if(silResolution)
                    eLoss1 = gRandom->Gaus(eLoss1, sigmaSil * TMath::Sqrt(eLoss1 / 5.5));
                T3AfterSil1 = T3AfterInterGas - eLoss1;
                ApplyNaN(eLoss1, specs.GetLayer(secondLayer.Data()).GetThresholds().at(silIndex1));
                isPunch = true;
            }
        }

        // Energy before the silicon
        double EBefSil0 {};
        if(isPunch && T3AfterSil1 == 0 && std::isfinite(eLoss1)) // if in f1
        {
            double EAfterSil0 {srim->EvalInitialEnergy("light", eLoss1, distance1)};
            EBefSil0 = eLoss0 + EAfterSil0;
        }
        else if(!isPunch && T3AfterSil0 == 0)
            EBefSil0 = eLoss0;
        else
            EBefSil0 = -1;

        ////////////////////////////////////////////////////////////////////////////////
        // 7->
        // we are ready to reconstruct Eex with all resolutions implemented
        //(d,light) is investigated gating on Esil1 = 0!

        // bool cutELoss0 {eLoss0Cuts[firstLayer].first <= eLoss0 && eLoss0 <= eLoss0Cuts[firstLayer].second}; // apply
        // pid cut
        bool cutEAfterSil0 {T3AfterSil0 == 0.}; // remove punchthrough

        // if(cutEAfterSil0 && std::isfinite(eLoss0) && cutELoss0) // if pid cuts applied
        if(cutEAfterSil0 && std::isfinite(EBefSil0)) // fill histograms
        {
            auto T3Rec {srim->EvalInitialEnergy("light", EBefSil0, distance0)};
            auto ExAfter {kin->ReconstructExcitationEnergy(T3Rec, theta3Lab)};
            auto thetaCM {kin->ReconstructTheta3CMFromLab(T3Rec, theta3Lab)};
            auto TBeamRec {kin->ReconstructBeamEnergyFromLabKinematics(T3Rec, theta3Lab)};
            // // Ecm from formula using masses
            auto ECMRec {(p2.GetAMU() / (p2.GetAMU() + p1.GetAMU())) * TBeamRec};
            auto TBeamDirRec {kin->ComputeEquivalentOtherT1(TBeamRec)};

            // fill histograms
            hECMRes->Fill(ECM, ECMRec);
            hTBeamDirRes->Fill(TBeamDir, TBeamDirRec);
            hKinSampled->Fill(theta3LabEff * TMath::RadToDeg(), T3Lab);
            hKinVertex->Fill(theta3Lab * TMath::RadToDeg(), T3Rec);
            hEexAfter->Fill(ExAfter, weight);
            hDistF0->Fill(distance0);
            hSilPID[layer0]->Fill(T3EnteringSil, (T3Lab - T3EnteringSil) / distance0);

            // Efficiency computation: passed histogram
            // WARNING: we must use the original thetaCM generated by the kinematic generator
            // otherwise we could be biasing the efficiency: a original thetaCM could be reconstructed shifted
            // (which makes sense after implementing resolutions) and hence contributing to another bin!!!
            // Besides, this could cause errors when making the division: passed counts > 0 / all counts == 0!!
            hThetaCM->Fill(thetaCMEff * TMath::RadToDeg());
            hTheta3Lab->Fill(theta3Lab * TMath::RadToDeg());
            hPhiLab->Fill(phi3Lab * TMath::RadToDeg());
            hEffAfter->Fill(thetaCMEff * TMath::RadToDeg(), ECM);
            hEffDirAfter->Fill(theta3LabDir,TBeamDir);
            hEffECNAfter->Fill(thetaCMEff * TMath::RadToDeg(),ECN);
            hECN->Fill(ECN);
            // ECN histograms from front and lateral, for comparison
            if(layer0 == "f0")
            {
                hECNFr->Fill(ECN);
                hSilSP[layer0]->Fill(silPoint0.Y(), silPoint0.Z());
                hTheta3CMfront->Fill(thetaCMEff);
                hTheta3Labfront->Fill(theta3Lab * TMath::RadToDeg());
                hRPz->Fill(vertex.Y(), vertex.Z());
            }
            else
            {
                hECNLat->Fill(ECN);
                hSilSP[layer0]->Fill(silPoint0.X(), silPoint0.Z());
                hTheta3CMside->Fill(thetaCMEff);
                hTheta3Labside->Fill(theta3Lab * TMath::RadToDeg());
            }

            // RP histogram
            hRP->Fill(vertex.X(), vertex.Y());
            // Beam energy at RP
            // hEBeam->Fill(theta3Lab*TMath::RadToDeg(),lightRange);

            // write to TTree
            Eex_tree = ExAfter;
            weight_tree = weight;
            theta3CM_tree = thetaCM * TMath::RadToDeg();
            EVertex_tree = T3Rec;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            rpx_tree = vertex.X();
            ecn_tree = ECN;
            hRPxSimu->Fill(vertex.X());
            outTree->Fill();
        }
        // hProj = hptoection->ProjectionX("", 50, 60);
    }
    std::cout << "\r" << std::string(100 / percentPrint, '|') << 100 << "%";
    std::cout.flush();
    std::cout << RESET << '\n';

    // Efficiencies as quotient of histograms in TEfficiency class
    auto* eff {new TEfficiency(*hThetaCM, *hThetaCMAll)};
    eff->SetNameTitle("eff", TString::Format("#theta_{CM} eff E_{x} = %.2f MeV", Ex));
    auto* effLab {new TEfficiency(*hThetaLab, *hThetaLabAll)};
    effLab->SetNameTitle("effLab", TString::Format("#theta_{Lab} eff E_{x} = %.2f MeV", Ex));
    auto* effPhi {new TEfficiency(*hPhiLab, *hPhiAll)};
    effPhi->SetNameTitle("effPhi", TString::Format("#phi_{Lab} eff E_{x} = %.2f MeV", Ex));
    // Manual computation of efficiencies
    auto* hEff2D {(TH2D*)hEffAfter->Clone("hEff2D")};
    hEff2D->Divide(hEffAll.get());
    auto* hEffDir2D {(TH2D*)hEffDirAfter->Clone("hEffDir2D")};
    hEffDir2D->Divide(hEffDirAll.get());
    auto* hEffECN2D {(TH2D*)hEffECNAfter->Clone("hEffECN2D")};
    hEffECN2D->Divide(hEffECNAll.get());

    // SAVING
    outFile->cd();
    outTree->Write();
    eff->Write();
    effLab->Write();
    effPhi->Write();
    for(auto& [_, h] : hSilSP)
        h->Write();
    hRP->Write("hRP");
    hEff2D->Write("hEff2D");
    hEffDir2D->Write("hEffDir2D");
    hEffECN2D->Write("hEffECN2D");
    hECMRes->Write("hECMRes");
    hTBeamDirRes->Write("hTBeamDirRes");
    outFile->Close();
    delete outFile;
    outFile = nullptr;


    // plotting
    if(standalone)
    {
        // draw theoretical kinematics
        ActPhysics::Kinematics theokin {p1, p2, p3, p4, T1 * p1.GetAMU(), Ex};
        auto* gtheo {theokin.GetKinematicLine3()};

        auto* c0 {new TCanvas("c0", "Canvas for inspection 0")};
        c0->DivideSquare(6);
        c0->cd(1);
        hThetaCM->DrawClone();
        c0->cd(2);
        hThetaCMAll->DrawClone();
        c0->cd(3);
        // hDistF0->DrawClone();
        hKinSampled->DrawClone("colz");
        gtheo->Draw("l");
        c0->cd(4);
        // hThetaCMAll->DrawClone();
        hRP->DrawClone("colz");
        c0->cd(5);
        hRPz->DrawClone("colz");
        c0->cd(6);

        auto* c1 {new TCanvas("c1", "PID")};
        c1->DivideSquare(4);
        auto* c2 {new TCanvas("c2", "SP")};
        c2->DivideSquare(4);
        int canvas {1};
        for(auto l : silLayers)
        {
            c1->cd(canvas);
            hSilPID[l]->Fit("fit", "Q");
            hSilPID[l]->DrawClone("colz");
            c2->cd(canvas);
            hSilSP[l]->DrawClone("colz");
            makeGrid(l, smAll);
            canvas++;
        }


        auto* c3 {new TCanvas {"c3", "2D Eff canvas"}};
        c3->DivideSquare(6);
        c3->cd(1);
        hEff2D->DrawClone("colz");
        c3->cd(2);
        hEffDir2D->DrawClone("colz");
        c3->cd(3);
        hECMRes->DrawClone("colz");
        c3->cd(4);
        hTBeamDirRes->DrawClone("colz");
        c3->cd(5);
        auto htot = (TH1*)hECN->DrawClone();
        hECNFr->SetLineColor(kMagenta);
        auto hfront = (TH1*)hECNFr->DrawClone("same");
        hECNLat->SetLineColor(kOrange);
        auto hlat = (TH1*)hECNLat->DrawClone("same");
        auto leg = new TLegend(0.6, 0.1, 0.9, 0.3);
        leg->AddEntry(htot, "All", "l");
        leg->AddEntry(hlat, "Lateral", "l");
        leg->AddEntry(hfront, "Front", "l");
        leg->Draw();

        c3->cd(6);
        hEffECN2D->DrawClone("colz");
    }

    delete srim;
    delete rand;

    timer.Stop();
    timer.Print();
}
