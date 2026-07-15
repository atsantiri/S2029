#include "ActColors.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActTPCParameters.h"
#include "ActSilSpecs.h"

#include "TCanvas.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLegend.h"
#include "ROOT/TThreadedObject.hxx"
#include "Math/Point3Dfwd.h"
#include "Rtypes.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <random>

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

double AngleWithNormal(const ROOT::Math::XYZVector &dir, const ROOT::Math::XYZVector &normal)
{
  auto dot{dir.Unit().Dot(normal.Unit())};
  return TMath::ACos(dot);
}

std::pair<ROOT::Math::XYZPoint, ROOT::Math::XYZPoint> SampleVertex(double meanZ, double sigmaZ, double meanY, double sigmaY, TH3D *h, double lengthX)
{

  // X is always common for both manners
  double Xstart{0};
  double Xrp{gRandom->Uniform() * lengthX};
  // Y depends completely on the method of calculation
  double Ystart{-1};
  double Yrp{-1};
  // Z of beam at entrance
  double Zstart{gRandom->Gaus(meanZ, sigmaZ)};
  double Zrp{-1};
  // Ystart in this case is sampled from the histogram itself!
  double thetaXY{};
  double thetaXZ{};

  // Two options depending on the existance of a emittance histogram or not
  if (h)
  {
    // Ystart in this case is sampled from the histogram itself!
    double thetaXY{};
    double thetaXZ{};
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
  ROOT::Math::XYZPoint start{Xstart, Ystart, Zstart};
  ROOT::Math::XYZPoint vertex{Xrp, Yrp, Zrp};
  return {std::move(start), std::move(vertex)};
}

void ApplyNaN(double &val, double thresh = 0, const std::string &comment = "stopped")
{
  if (val <= thresh)
    val = std::nan(comment.c_str());
}

void makeGrid(std::string layer, std::unordered_map<std::string, ActPhysics::SilMatrix *> smAll)
{
  auto &cuts = smAll[layer]->GetGraphs();
  for (const auto &[id, cut] : cuts)
  {
    if (cut)
      cut->Draw("same");
  }
}

void Simulation_S2029(const std::string &beam, const std::string &target,
                      const std::string &light, const std::string &heavy,
                      double T1, double Ex, int pressure, bool standalone)
{
  // set batch mode if not an independent function
  if (!standalone)
    gROOT->SetBatch(true);

  double angle_min{130.};
  double angle_max{140.};
  double qval{3.923}; // qvalue of 17F+p

  TRandom random;

  // SIGMAS
  const double sigmaSil{0.05}; // AT: To change, get sigma from calibrations
  const double sigmaPercentBeam{0};
  const double sigmaAngleLight{0.95 / 2.355};
  // Parameters of beam in mm
  // Beam has to be manually placed in the simulation
  // Centered in Z and Y with a width of 4 mm
  // Center in Z
  // AT: note that in simue756 zVertexMean is coming from silicon matrices
  const double zVertexMean{128. + 18.}; // beam not centered in chamber, upwards
                                        //    double zVertexMean {silCentre + beamOffset}; // in Miguel's

  const double zVertexSigma{4}; // AT, check why diff from e796
  // Center in Y - miguel's simu doesn't have Yvertex?
  const double yVertexMean{128.};
  const double yVertexSigma{4};
  // const double zVertexMean {83.59};
  // const double zVertexSigma {3.79};

  // THRESHOLDS FOR SILICONS -- AT double check from calibration files
  const double thresholdSi0{1.};
  const double thresholdSi1{1.};

  // number of iterations
  const int iterations{static_cast<int>(1e6)};

  // ACTIVATE STRAGGLING OR NOT
  bool stragglingInGas{false};
  bool stragglingInSil{false};
  bool silResolution{false};
  bool thetaResolution{false};

  // TPC basic parameters
  ActRoot::TPCParameters fTPC{"Actar"};

  // Silicon specs
  ActPhysics::SilSpecs specs;
  specs.ReadFile("../configs/silspecs.conf");

  // Silicon EFFECTIVE matrix
  double silCentre{};
  double beamOffset{}; // determined from emittance calculations
  // std::string layer = "f0"; // for now
  // ActPhysics::SilMatrix *sm{specs.GetLayer(layer).GetSilMatrix()->Clone()};
  std::vector<std::string> silLayers{"f0", "l0", "r0"};
  std::unordered_map<std::string, ActPhysics::SilMatrix *> smAll;
  for (const auto &l : silLayers)
    smAll[l] = specs.GetLayer(l).GetSilMatrix()->Clone();

  /////////////////////////////////////////////////////////////////////////////
  // need to move silicons around, like Ivan's lines 304 - 334 of s2384/Simulation/do_simu.cxx
  /////////////////////////////////////////////////////////////////////////////
  specs.DrawGeo();

  // Init particles
  ActPhysics::Particle p1{beam};
  ActPhysics::Particle p2{target};
  ActPhysics::Particle p3{light};
  ActPhysics::Particle p4{heavy};

  // Init kinematics generator
  ActSim::KinematicGenerator kingen{p1, p2, p3, p4, 0, 0};
  kingen.Print();
  ActPhysics::Kinematics *reckin = kingen.GetBinaryKinematics();

  // Get threshold
  auto T1Thresh{ActPhysics::Kinematics(p1, p2, p3, p4, -1, Ex).GetT1Thresh()};

  // Histograms
  // To compute a fine-grain efficiency, we require at least a binning width of 0.25 degrees!
  auto hThetaCM{HistConfig::ThetaCM.GetHistogram()};
  auto hThetaCMAll{HistConfig::ChangeTitle(HistConfig::ThetaCM, "ThetaCM all", "All").GetHistogram()};
  auto hDistL0{HistConfig::ChangeTitle(HistConfig::TL, "Distance to L0").GetHistogram()};
  auto hKinVertex{HistConfig::ChangeTitle(HistConfig::Kin, "Kinematics at vertex").GetHistogram()};
  auto hEexBefore{HistConfig::ChangeTitle(HistConfig::Ex, "Ex before resolutions", "Bef").GetHistogram()};
  auto hEexAfter{HistConfig::ChangeTitle(HistConfig::Ex, "Ex after resolutions", "After").GetHistogram()};
  auto hSPTheta{std::make_unique<TProfile2D>("hSPTheta", "SP vs #theta_{CM};Y [mm];Z [mm];#theta_{CM} [#circ]", 75, 0, 300, 75, 0, 300)};
  auto hRP{HistConfig::RP.GetHistogram()};
  auto hEcn{HistConfig::ECN.GetHistogram()};
  auto hEcnFr{HistConfig::ECN.GetHistogram()};  // front Sil
  auto hEcnLat{HistConfig::ECN.GetHistogram()}; // lateral sil
  auto *hRPxSimu{HistConfig::RP.GetHistogram()->ProjectionX("hRPxSimu")};
  auto *hRPxRange{new TH2D{"hRPxRange", "Light range vs dist to sil wall;Dist to wall [mm];Light range [mm]", 300, 0, 500, 300, 0, 500}};
  auto *hEBeam{new TH2D{"hEBeam", "Range angle", 300, 0, 180, 300, 0, 500}};
  auto *hEx2{new TH2D{"hEx2", "#theta_{Lab} vs E_{CN};E_{CN} [MeV];#theta_{CM} [#circ]", 200, 0, 9, 200, 90, 180}};
  // auto *hnorm{new TH2D{"hnorm", "geometric efficiency;E*_{^{17}F,elastic} [MeV];#theta_{CM} [#circ]", 90, 4, 9, 90, 90, 180}};
  auto hnorm{HistConfig::ChangeTitle(HistConfig::EcmThetaCM, "geometric efficiency").GetHistogram()};
  auto hKin{HistConfig::Kin.GetHistogram()};
  // auto *hnorm_phi{new TH2D{"hnorm_phi", "geometric efficiency;E*_{^{17}F} [MeV];#Phi_{Lab} [#circ]", 90, 0, 9, 90, 0, 200}};
  auto hnorm_phi{HistConfig::ChangeTitle(HistConfig::EcmPhiLab, "geometric efficiency").GetHistogram()};
  auto hECM{HistConfig::ECM.GetHistogram()};
  auto *hprojection{new TH1D{"hprojection", "Ex;E_{x} [MeV]; Counts", 80, 11, 13}};
  auto *htrack{new TH2D{"htrack", "trackslength;#theta_{heavy lab} [#circ];track length [mm]", 90, 0, 100, 90, 0, 250}};
  auto hEcnThetaCM{HistConfig::EcnThetaCM.GetHistogram()};
  auto *hTestEcm{new TH2D{"hTestEcm", "E_{CM};E_{CM} ;E_{CM}Recon", 100, 0, 5, 100, 0, 5}};
  auto hTheta3CM{HistConfig::ThetaCM.GetHistogram()};
  auto hTheta3CMside{HistConfig::ThetaCM.GetHistogram()};
  auto hTheta3CMfront{HistConfig::ThetaCM.GetHistogram()};
  auto hTheta3Lab{HistConfig::ThetaCM.GetHistogram()};
  hTheta3Lab->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
  auto hTheta3Labside{HistConfig::ThetaCM.GetHistogram()};
  hTheta3Labside->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
  auto hTheta3Labfront{HistConfig::ThetaCM.GetHistogram()};
  hTheta3Labfront->SetTitle("Theta3Lab;#theta_{Lab} [#circ];Counts");
  auto hThetaLabAll{HistConfig::ThetaCM.GetHistogram()};
  hThetaLabAll->SetTitle("Theta3Lab all;#theta_{Lab} [#circ];Counts");
  auto hPhiAll{HistConfig::PhiCM.GetHistogram()};
  hPhiAll->SetTitle("Phi3CM all;#phi_{CM} [#circ];Counts");

  std::map<std::string, ROOT::TThreadedObject<TH2D>> hSilPID, hSilSP;
  auto hSP{HistConfig::SP.GetHistogram()};
  auto *hpid{new TH2D{"hPID", "PID;E[MeV];Eloss [MeV]", 1000, 0, 10, 1000, 0, 0.01}};
  for (const auto &layer : {"f0", "l0", "r0"})
  {
    hSilPID.emplace(layer, *hpid);
    hSilPID[layer]->SetTitle(TString::Format("%s", layer));
    hSilSP.emplace(layer, *hSP);
    hSilSP[layer]->SetTitle(TString::Format("%s", layer));
  }
  auto hSPCut{(TH2D *)hSP->Clone("hSPCut")}; // 3D object to examine L1 events?

  // Load SRIM tables
  // The name of the file sets particle + medium
  auto *srim{new ActPhysics::SRIM()};
  srim->ReadTable("light", TString::Format("SRIM/%s_H2-iC4H10_95-5_%dmbar.txt", light.c_str(), pressure).Data());
  srim->ReadTable("beam", TString::Format("SRIM/%s_H2-iC4H10_95-5_%dmbar.txt", beam.c_str(), pressure).Data());
  srim->ReadTable("lightInSil", TString::Format("SRIM/%s_silicon.txt", light.c_str()).Data());

  // Random generator
  auto *rand{new TRandom3()};
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
  auto *outFile{new TFile(TString::Format("../PostAnalysis/Input/Simu_%s_%dmbar.root", light.c_str(), pressure), "recreate")};
  auto *outTree{new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};

  double theta3CM_tree{};
  outTree->Branch("theta3CM", &theta3CM_tree);
  double Eex_tree{};
  outTree->Branch("Eex", &Eex_tree);
  double weight_tree{};
  outTree->Branch("weight", &weight_tree);
  double EVertex_tree{};
  outTree->Branch("EVertex", &EVertex_tree);
  double theta3Lab_tree{};
  outTree->Branch("theta3Lab", &theta3Lab_tree);
  double rpx_tree{};
  outTree->Branch("RPx", &rpx_tree);
  double ecn_tree{};
  outTree->Branch("ECN", &ecn_tree);

  //---- SIMULATION STARTS HERE
  ROOT::EnableImplicitMT();

  // timer
  TStopwatch timer{};
  timer.Start();

  // print fancy info
  std::cout << BOLDMAGENTA << "Running for Ex = " << Ex << " MeV" << RESET << '\n';
  std::cout << BOLDMAGENTA << "Output File: " << outFile->GetName() << RESET << '\n';
  std::cout << BOLDGREEN;
  const int percentPrint{5};
  int step{iterations / (100 / percentPrint)};
  int nextPrint{step};
  int percent{};
  double maxBinContent = -1;
  for (long int reaction = 0; reaction < iterations; reaction++)
  {
    // Print progress
    if (reaction >= nextPrint)
    {
      percent = 100 * (reaction + 1) / iterations;
      int nchar{percent / percentPrint};
      std::cout << "\r" << std::string((int)(percent / percentPrint), '|')
                << percent << "%";
      std::cout.flush();
      nextPrint += step;
    }
    // 1-> Sample vertex
    auto [start, vertex]{SampleVertex(zVertexMean, zVertexSigma, yVertexMean, yVertexSigma, nullptr, fTPC.X())};

    // 2-> Beam energy according to its sigma
    auto TBeam{runner.RandomizeBeamEnergy(T1 * p1.GetAMU(), sigmaPercentBeam * T1 * p1.GetAMU())}; // T1 in Mev / u * mass of beam in u = total kinetic energy

    // Slow down it according to vertex position
    auto distToVertex{(vertex - start).R()};
    // std::cout << "distToVertex = " << distToVertex << "\n";

    // If distance is too small, protect from spline faults
    srim->SetUseSpline(false);
    TBeam = srim->Slow("beam", TBeam, distToVertex);
    srim->SetUseSpline(true);

    // runner energy functions return std::nan when the particle is stopped in
    // the gas! if nan (aka stopped in gas, continue) if not stopped but beam
    // energy below kinematic threshold, continue
    if (std::isnan(TBeam) || TBeam < T1Thresh)
      continue;

    ////////////////////////////////////////////////////////////////////////////////
    // 3-> Run kinematics!
    kingen.SetBeamAndExEnergies(TBeam, Ex);
    reckin->SetBeamEnergyAndEx(TBeam, Ex);
    double weight{kingen.Generate()}; // need to understand what this is

    // Light
    auto *PLight{kingen.GetLorentzVector(0)};
    auto theta3Lab{PLight->Theta()};
    auto phi3Lab{PLight->Phi()};
    auto T3Lab{PLight->Energy() - p3.GetMass()};

    // Heavy
    auto *PHeavy{kingen.GetLorentzVector(1)};
    auto T4Lab{PHeavy->Energy() - p4.GetMass()};

    // Save without resolution
    // Simulated thetaCM, to be used for efficiency computation by CM interval and with our set reference direction
    double theta3CMEff{reckin->ReconstructTheta3CMFromLab(T3Lab, theta3Lab)};
    double phi3CM = phi3Lab;
    auto EexEff{reckin->ReconstructExcitationEnergy(T3Lab, theta3Lab)};
    double theta3LabEff{theta3Lab}; // before implementing resolution in angle
    auto Ecm{(p2.GetAMU() / (p2.GetAMU() + p1.GetAMU())) * TBeam};

    // Fill Histograms
    hThetaCMAll->Fill(theta3CMEff * TMath::RadToDeg());
    hThetaLabAll->Fill(theta3LabEff * TMath::RadToDeg());
    hPhiAll->Fill(phi3Lab * TMath::RadToDeg());

    ////////////////////////////////////////////////////////////////////////////////
    // 4-> Include thetaLab resolution to compute thetaCM and Ex afterwards
    if (thetaResolution) // resolution in
      theta3Lab = gRandom->Gaus(theta3Lab, sigmaAngleLight * TMath::DegToRad());

    // Eval range of light particle
    auto lightRange{srim->EvalDirect("light", T3Lab)};
    // std::cout<<"lightRange = "<<lightRange<<std::endl;
    // auto heavyRange{srim->EvalDirect("heavy", T4Lab)};

    ////////////////////////////////////////////////////////////////////////////////
    // 5-> Propagate track from vertex to silicon wall using SilSpecs class
    // And using the angle with the uncertainty already in
    ROOT::Math::XYZVector dirBeamFrame{TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab), TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};

    // Declare beam direction
    auto beamDir{(vertex - start).Unit()};
    // Rotate to world = geometry frame
    auto dirWorldFrame{runner.RotateToWorldFrame(dirBeamFrame, beamDir)};

    // Threshold L1, particles that stop in actar. Check before doing the continues
    double rangeInGas{srim->EvalRange("light", T3Lab)};
    ROOT::Math::XYZPoint finalPointGas{vertex + rangeInGas * dirWorldFrame.Unit()};
    if (0 <= finalPointGas.X() && finalPointGas.X() <= fTPC.X() && 0 <= finalPointGas.Y() && finalPointGas.Y() <= fTPC.Y() && 0 <= finalPointGas.Z() && finalPointGas.Z() <= fTPC.Z())
    {
      hSPCut->Fill(finalPointGas.Y(), finalPointGas.Z());
      // do something here????
    }

    // auto [silIndex0, silPoint0InMM]{specs.FindSPInLayer(layer, vertex, dirWorldFrame)};
    int silIndex0 = -1;
    ROOT::Math::XYZPoint silPoint0;
    std::string layer0;
    for (auto layer : silLayers)
    {
      std::tie(silIndex0, silPoint0) = specs.FindSPInLayer(layer, vertex, dirWorldFrame);
      if (silIndex0 != -1)
      {
        layer0 = layer;
        break;
      }
    }
    // Define SP distance
    auto distance0{(vertex - silPoint0).R()};
    // std::cout << "lightRange: " << lightRange << " distance0: " << distance0 << std::endl;
    auto T3EnteringSil{srim->Slow("light", T3Lab, distance0)};
    // std::cout<<" SilPoint "<<silPoint0InMM<<", energy entering Sil: "<<T3EnteringSil<<std::endl;
    ApplyNaN(T3EnteringSil);
    // nan if stopped in gas
    if (!std::isfinite(T3EnteringSil))
      continue;

    // skip tracks that doesn't reach silicons or are in silicon index cut
    if (silIndex0 == -1)
    {
      continue;
    }

    // SILICON0
    ROOT::Math::XYZVector normal{specs.GetLayer(layer0).GetNormal()};
    auto angleNormal0{AngleWithNormal(dirWorldFrame, normal)};
    auto T3AfterSil0{srim->SlowWithStraggling("lightInSil", T3EnteringSil, specs.GetLayer(layer0).GetUnit().GetThickness(), 0)}; ////////////////////////////////////////////////////////////////////////////////////////////make angleWithNormal
    auto eLoss0{T3EnteringSil - T3AfterSil0};

    // Apply resolution
    if (T3AfterSil0 != 0)
    {
      eLoss0 = gRandom->Gaus(eLoss0, sigmaSil * TMath::Sqrt(eLoss0 / 5.5));
      T3AfterSil0 = T3EnteringSil - eLoss0;
    }
    // ApplyNaN(eLoss0, thresholdSi0, "thresh");
    ApplyNaN(eLoss0, specs.GetLayer(layer0).GetThresholds().at(silIndex0));

    // nan if bellow threshold
    if (!std::isfinite(eLoss0))
      continue;

    ////////////////////////////////////////////////////////////////////////////////
    // 6-> Same but to silicon layer 1 if exists

    ////////////////////////////////////////////////////////////////////////////////
    // 7->
    // we are ready to reconstruct Eex with all resolutions implemented
    //(d,light) is investigated gating on Esil1 = 0!
    bool cutEAfterSil0{T3AfterSil0 == 0.};      // remove punchthrough
    if (cutEAfterSil0 && std::isfinite(eLoss0)) // fill histograms
    {
      // auto T3Recon{runner.EnergyBeforeGas(eLoss0, distance0, "light")};
      auto T3Rec{srim->EvalInitialEnergy("light", eLoss0, distance0)};
      auto EexRec{reckin->ReconstructExcitationEnergy(T3Rec, theta3Lab)};
      auto theta3CMRec{reckin->ReconstructTheta3CMFromLab(T3Rec, theta3Lab)};
      auto TBeamRec{reckin->ReconstructBeamEnergyFromLabKinematics(T3Rec, theta3Lab)};
      // // Ecm from formula using masses
      auto EcmRecon{(p2.GetAMU() / (p2.GetAMU() + p1.GetAMU())) * TBeamRec};
      double EcnRec = Ecm + qval;
      // double theta_gauss = random.Gaus(theta3Lab, 0.1);

      auto Ex_calc = T3Lab * ((2 * TMath::Cos(theta3Lab) * TMath::Sqrt((TBeam * p3.GetAMU()) / (T3Lab * p4.GetAMU()))) - ((p3.GetAMU() / p4.GetAMU()) + 1));

      // fill histograms
      hTheta3CM->Fill(theta3CMEff * TMath::RadToDeg());
      hThetaCM->Fill(theta3CMRec * TMath::RadToDeg());
      TF1 *fit = new TF1("fit", "[0]*(x^[1])", 10, -0.5);
      hTheta3Lab->Fill(theta3Lab * TMath::RadToDeg());
      hSilPID[layer0]->Fill(T3EnteringSil, (T3Lab - T3EnteringSil) / distance0);
      // hThetaCM->Fill(theta3CMBefore * TMath::RadToDeg());
      hEexBefore->Fill(EexEff, weight); // with the weight from each TGenPhaseSpace::Generate()
      hDistL0->Fill(distance0);
      hKin->Fill(theta3Lab * TMath::RadToDeg(), T3EnteringSil);
      hKinVertex->Fill(theta3Lab * TMath::RadToDeg(), Ecm); // T3Recon);
      hEcnThetaCM->Fill(EcnRec, theta3CMRec * TMath::RadToDeg());
      hnorm->Fill(EcmRecon, theta3CMRec * TMath::RadToDeg());
      hnorm_phi->Fill(EcmRecon, phi3Lab * TMath::RadToDeg());
      hEcn->Fill(EcnRec);
      if (layer0 == "f0")
      {
        hEcnFr->Fill(EcnRec);
        hSilSP[layer0]->Fill(silPoint0.Y(), silPoint0.Z());
        hTheta3CMfront->Fill(theta3CMEff);
        hTheta3Labfront->Fill(theta3Lab * TMath::RadToDeg());
      }
      else
      {
        hEcnLat->Fill(EcnRec);
        hSilSP[layer0]->Fill(silPoint0.X(), silPoint0.Z());
        hTheta3CMside->Fill(theta3CMEff);
        hTheta3Labside->Fill(theta3Lab * TMath::RadToDeg());
      }
      hTestEcm->Fill(Ecm, EcmRecon);

      hEexAfter->Fill(EexRec, weight);
      // hSP->Fill(silPoint0.Y(), silPoint0.Z());

      // Fill histogram of SP with thetaCM as weight
      hSPTheta->Fill(silPoint0.Y(), silPoint0.Z(), theta3CMRec * TMath::RadToDeg());
      // RP histogram
      hRP->Fill(vertex.X(), vertex.Y());
      // Beam energy at RP
      // hEBeam->Fill(theta3Lab*TMath::RadToDeg(),lightRange);
      // CM energy
      hECM->Fill(Ecm);

      // write to TTree
      Eex_tree = EexRec;
      weight_tree = weight;
      theta3CM_tree = theta3CMRec * TMath::RadToDeg();
      EVertex_tree = T3Rec;
      theta3Lab_tree = theta3Lab * TMath::RadToDeg();
      rpx_tree = vertex.X();
      ecn_tree = EcnRec;
      hRPxSimu->Fill(vertex.X());
      outTree->Fill();
    }
    // hProj = hptoection->ProjectionX("", 50, 60);
  }

  // Obtenir le nombre de bins
  int nBinsX = hnorm->GetNbinsX();
  int nBinsY = hnorm->GetNbinsY();

  for (int i = 1; i <= nBinsX; i++)
  {
    for (int j = 1; j <= nBinsY; j++)
    {

      double n_observe = hnorm->GetBinContent(i, j);
      // std::cout << "n_data" << n_observe << "\n";

      // double n_theorique = hnorm_uniform->GetBinContent(i, j);
      double n_theorique = (iterations / 8100.) * 3.;
      // std::cout << "n_norm" << n_theorique << "\n";
      // std::cout << "\n";

      // Normalise le bin
      double normalise = 0;
      if (n_theorique > 0)
      {
        normalise = n_observe / n_theorique;
      }

      // Mets à jour le contenu du bin normalisé
      hnorm->SetBinContent(i, j, normalise);
    }
  }

  // Compute efficiency side, front, and total
  auto *effCM{new TEfficiency{*hTheta3CM, *hThetaCMAll}};
  effCM->SetNameTitle("effCM", " #epsilon_{TOT} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
  auto *effLab{new TEfficiency{*hTheta3Lab, *hThetaLabAll}};
  effLab->SetNameTitle("effLab", "#epsilon_{TOT} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

  auto *effCMside{new TEfficiency{*hTheta3CMside, *hThetaCMAll}};
  effCMside->SetNameTitle("effCMside", " #epsilon_{side} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
  auto *effLabside{new TEfficiency{*hTheta3Labside, *hThetaLabAll}};
  effLabside->SetNameTitle("effLabside", "#epsilon_{side} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

  auto *effCMfront{new TEfficiency{*hTheta3CMfront, *hThetaCMAll}};
  effCMfront->SetNameTitle("effCMfront", " #epsilon_{front} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
  auto *effLabfront{new TEfficiency{*hTheta3Labfront, *hThetaLabAll}};
  effLabfront->SetNameTitle("effLabfront", "#epsilon_{front} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");

  // SAVING
  outFile->cd();
  hnorm->Write();
  hnorm_phi->Write();
  hprojection->Write();
  outTree->Write();
  effCM->Write();
  effLab->Write();
  effCMside->Write();
  effLabside->Write();
  effCMfront->Write();
  effLabfront->Write();
  outFile->Close();
  delete outFile;
  outFile = nullptr;

  // plotting
  if (standalone)
  {
    // draw theoretical kinematics
    ActPhysics::Kinematics theokin{p1, p2, p3, p4, T1 * p1.GetAMU(), Ex};
    auto *gtheo{theokin.GetKinematicLine3()};

    set_plot_style();
    auto *c0{new TCanvas("c0", "PID")};
    c0->DivideSquare(4);
    auto *c2{new TCanvas("c2", "SP")};
    c2->DivideSquare(4);
    int canvas{1};
    for (auto l : silLayers)
    {
      c0->cd(canvas);
      hSilPID[l]->Fit("fit");
      hSilPID[l]->DrawClone("colz");
      c2->cd(canvas);
      hSilSP[l]->DrawClone("colz");
      makeGrid(l, smAll);
      canvas++;
    }
    auto *c1{new TCanvas("cAfter", "Canvas for inspection 1")};
    c1->DivideSquare(6);
    c1->cd(1);
    hnorm->DrawClone("colz");
    c1->cd(2);
    // effCM->SetTitle("SP for particles not reaching sils");
    hThetaCMAll->DrawClone();
    c1->cd(5);
    hTheta3CM->DrawClone();
    c1->cd(3);
    auto htot = (TH1 *)hEcn->DrawClone();
    hEcnFr->SetLineColor(kMagenta);
    auto hfront = (TH1 *)hEcnFr->DrawClone("same");
    hEcnLat->SetLineColor(kOrange);
    auto hlat = (TH1 *)hEcnLat->DrawClone("same");
    auto leg = new TLegend(0.6, 0.1, 0.9, 0.3);
    leg->AddEntry(htot, "All", "l");
    leg->AddEntry(hlat, "Lateral", "l");
    leg->AddEntry(hfront, "Front", "l");
    leg->Draw();

    c1->cd(4);
    hEcnThetaCM->DrawClone("col");
  }

  delete srim;
  delete rand;

  timer.Stop();
  timer.Print();
}
