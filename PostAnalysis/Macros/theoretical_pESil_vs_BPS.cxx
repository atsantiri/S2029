// compute theoretical line for pESil_vs_BSP.cxx
#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActRunner.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"

#include "TCanvas.h"
#include "TMath.h"

// function to propagate proton track from vertex to Si wall and return ESil - based on Simulation_E796.cpp
// note that I have full paths instead of relative paths, as this macro can be executed from different scripts
double calc_ESil(double RPx, double theta3Lab, double T3Lab)
{
    // Silicon specs
    ActPhysics::SilSpecs specs;
    specs.ReadFile("../configs/silspecs.conf");
    // specs.ReadFile("/home/artemis/ACTAR/S2029/configs/silspecs.conf");

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("gas", "../Simulation/SRIM/1H_H2-iC4H10_95-5_760mbar.txt");
    // srim->ReadTable("gas", "/home/artemis/ACTAR/S2029//Simulation/SRIM/1H_H2-iC4H10_95-5_760mbar.txt");
    srim->ReadTable("silicon", "../Simulation/SRIM/1H_silicon.txt");
    // srim->ReadTable("silicon", "/home/artemis/ACTAR/S2029//Simulation/SRIM/1H_silicon.txt");

    // Runner: contains utility functions to execute multiple actions
    ActSim::Runner runner(nullptr, nullptr, 0, 0.);

    // specify track frame relative to the beam frame
    double phi {0};                               // make life easy for now, will revisit later if needed
    ROOT::Math::XYZVector beamDir(1.0, 0.0, 0.0); // assume beam goes along the x axis for simplicity
    ROOT::Math::XYZPoint vertex(RPx, 122, 135);   // assume beam is centered on Z, Y is an approximation from what I saw in Macros/Beam/beamEmittance.cxx

    ROOT::Math::XYZVector dirBeamFrame {TMath::Cos(theta3Lab * TMath::DegToRad()),
                                        TMath::Sin(theta3Lab * TMath::DegToRad()) * TMath::Sin(phi),
                                        TMath::Sin(theta3Lab * TMath::DegToRad()) * TMath::Cos(phi)};
    // Rotate to world = geometry frame
    auto dirWorldFrame {runner.RotateToWorldFrame(dirBeamFrame, beamDir)};

    // Find which silicon will it hit, if any.
    std::string siLayer {};
    int silIdx {-1};
    ROOT::Math::XYZPoint hitPoint {};
    bool flag {false};
    for(const auto& layer : {"f0", "l0", "r0"})
    {
        auto [silIndex0, silPoint0InMM] {specs.FindSPInLayer(layer, vertex, dirWorldFrame)};
        if(silIndex0 != -1)
        {
            if(flag)
                throw std::runtime_error("I'm confused! Particle hit multiple silicons...");
            silIdx = silIndex0;
            hitPoint = silPoint0InMM;
            siLayer = layer;
            flag = true;
        }
    }

    if(silIdx == -1)
        return 0.;

    // Propagate particle to impact point and find deposited energy
    auto distance {(vertex - hitPoint).R()};
    auto T3EnteringSil {srim->Slow("gas", T3Lab, distance)};

    ROOT::Math::XYZVector silDir = (siLayer == "f0") ? ROOT::Math::XYZVector(1, 0, 0) : ROOT::Math::XYZVector(0, 1, 0);
    auto angleNormal {TMath::ACos(dirWorldFrame.Unit().Dot(silDir))};

    auto T3AfterSil {
        srim->Slow("silicon", T3EnteringSil, specs.GetLayer(siLayer).GetUnit().GetThickness(), angleNormal)};
    auto eLoss {T3EnteringSil - T3AfterSil};

    return eLoss;
}


// create theoretical line for given lab angle and Eex
TGraph* calcTheo_pESil_vs_BSP(double theta3Lab, double Eex, EColor color, int style)
{
    TGraph* ret = new TGraph();
    ret->SetTitle(TString::Format("Ex=%.02f MeV", Eex));
    ret->SetLineWidth(2);
    ret->SetLineColor(color);
    ret->SetLineStyle(style);

    // Init SRIM
    auto* srim {new ActPhysics::SRIM};
    // Correct SRIM names
    std::string light = "1H";
    std::string beam = "17F";
    srim->ReadTable(beam, TString::Format("../Simulation/SRIM/%s_H2-iC4H10_95-5_760mbar.txt", beam.c_str()).Data());
    // srim->ReadTable(beam, TString::Format("/home/artemis/ACTAR/S2029/Simulation/SRIM/%s_H2-iC4H10_95-5_760mbar.txt", beam.c_str()).Data());

    // Init particles
    ActPhysics::Particle pb {beam};
    ActPhysics::Particle pt {light};
    ActPhysics::Particle pl {light};

    // Initial energy of beam at pad plane entrance
    double EBeamIni {3.84}; // MeV/u

    // calculate reaction threshold energy
    ActPhysics::Kinematics kin {pb, pt, pt, EBeamIni * pb.GetAMU(), Eex};
    auto beamThreshold {kin.GetT1Thresh()};

    // loop through possible RPx distances
    double xmin {0.};
    double xmax {srim->EvalRange(beam, EBeamIni * pb.GetAMU())};
    int nsteps {200};

    // std::cout << "x \tEbeam \tthCM \tth4L \tT4Lab \trT4 \tT3Lab \tBSP" << std::endl;
    for(int i = 0; i < nsteps; i++)
    {
        double x = xmin + i * (xmax - xmin) / nsteps;
        // Propagate beam to RPx
        double Ebeam {srim->Slow(beam, EBeamIni * pb.GetAMU(), x)};

        // If reaction possible calculate kinematics
        if(Ebeam > beamThreshold)
        {
            ActPhysics::Kinematics kin {pb, pt, pt, Ebeam, Eex};

            double thetaCM {kin.ReconstructTheta3CMFromLab(Ebeam, theta3Lab * TMath::DegToRad()) * TMath::RadToDeg()};
            kin.ComputeRecoilKinematics(thetaCM * TMath::DegToRad(), 0.);

            double theta4Lab {kin.GetTheta4Lab() * TMath::RadToDeg()};
            double T3Lab {kin.GetT3Lab()};
            double T4Lab {kin.GetT4Lab()};
            double rangeT4 {srim->EvalRange(beam, T4Lab)};

            // Beam-like stopping point is RPx + TLheavy
            double BSP {x + rangeT4};

            double ESil {calc_ESil(x, theta3Lab, T3Lab)};
            if (ESil==0)
                continue;

            // std::cout << std::setprecision(3) << x << "\t" << Ebeam << "\t" << T4Lab / pb.GetAMU() << "\t" << rangeT4
            //           << "\t" << T3Lab << "\t" << BSP << "\t" << ESil << std::endl;


            ret->SetPoint(ret->GetN(), BSP, ESil);
        }
    }


    return ret;
}

void theoretical_pESil_vs_BPS()
{
    auto* c0 {new TCanvas("c0", "ESil vs BSP")};
    auto gr {calcTheo_pESil_vs_BSP(17.5, 0., kRed, 1)};
    gr->Draw();
}