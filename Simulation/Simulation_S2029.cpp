// Based Simulation_E796.cpp

#include "ActColors.h"
// #include "ActCrossSection.h"
// #include "ActDecayGenerator.h"
// #include "ActKinematicGenerator.h"
// #include "ActKinematics.h"
// #include "ActLine.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
// #include "ActUtils.h"

#include <TCanvas.h>
#include <TLegend.h>


#include "../PostAnalysis/HistConfig.h"

double dEthroughIC(const ActPhysics::Particle& particle, double energy, bool applyStraggling, ActPhysics::SRIM* srim) {
    if (!applyStraggling) return energy;

    // Choose the right SRIM table key based on particle name
    TString tableKey1, tableKey2;
    if (particle.GetName() == "17F" || particle.GetName() == "20Mg") {
        tableKey1 = "beamInMylar";
        tableKey2 = "beamInCF4";
    } else if (particle.GetName() == "17O" || particle.GetName() == "20Na") {
        tableKey1 = "ContaminantInMylar";
        tableKey2 = "ContaminantInCF4";
    } else {
        std::cerr << "Unknown particle in dEthroughIC: " << particle.GetName() << '\n';
        return energy; 
    }
    
    // 11.5 um mylar
    auto eneAfterICmylar {srim->SlowWithStraggling(tableKey1.Data(), energy, 11.5e-3)}; // 11.5 um mylar in IC 
    // 36 mm CF4 gas
    auto eneAfterIC {srim->SlowWithStraggling(tableKey2.Data(), eneAfterICmylar, 36)};
    
    return eneAfterIC;
}

void Simulation_S2029(const std::string& beam = "17F", double T1 = 4.5, double Ex = 0){

    // Set number of iterations
    const int iterations {static_cast<int>(1e6)};
    gRandom->SetSeed();
    // Runner: contains utility functions to execute multiple actions
    ActSim::Runner runner(nullptr, nullptr, gRandom, 0);


    // Set parameters to include
    bool stragglingInIC {true};   // If true beam will be propagated through IC
    double pressureIC {70}; // Pressure in mbar    
    bool stragglingInCFA {true};  // If true beam will be propagated through CFA
    double pressureCFA {6}; // Pressure in mbar    
    bool trackContamination {true}; // If true the 17O contaminant will be propagated along with the 17F primary beam

    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silspecs.conf");
    sils->Print();

    // Resolutions
    const double sigmaSil {0.060 / 2.355};
    const double sigmaPercentBeam {0.004};
    const double sigmaAngleLight {0.95 / 2.355};

    // Kinematics
    ActPhysics::Particle p1 {beam};
    ActPhysics::Particle p1c {"17O"}; // Contaminant 
    // ActPhysics::Particle p2 {target};
    // ActPhysics::Particle p3 {arglight};
    // // Automatically compute 4th particle
    // ActPhysics::Kinematics kaux {p1, p2, p3};
    // ActPhysics::Particle p4 {kaux.GetParticle(4)};
    // // Binary kinematics generator
    // ActSim::KinematicGenerator kingen {
    //     p1, p2, p3, p4, (protonPS > 0 ? protonPS : 0), (neutronPS > 0 ? neutronPS : (pdphase || dtphase ? 1 : 0))};
    // kingen.Print();

    // Load SRIM Tables
    auto* srim {new ActPhysics::SRIM()};
    TString srim_IC_mylar = TString::Format("../Calibrations/SRIM/%s_mylar.txt", beam.c_str());
    if (std::ifstream(srim_IC_mylar.Data())) {
        srim->ReadTable("beamInMylar", srim_IC_mylar.Data());
    } else {
        std::cerr << srim_IC_mylar <<" not found!\n";
        return;
    }
    TString srim_IC_gas = TString::Format("../Calibrations/SRIM/%s_CF4_%0.fmbar.txt", beam.c_str(), pressureIC);
    if (std::ifstream(srim_IC_gas.Data())) {
        srim->ReadTable("beamInCF4", srim_IC_gas.Data());
    } else {
        std::cerr << srim_IC_gas <<" not found!\n";
        return;
    }
    TString srim_p1c_mylar = TString::Format("../Calibrations/SRIM/%s_mylar.txt", p1c.GetName().c_str());
    if (std::ifstream(srim_p1c_mylar.Data())) {
        srim->ReadTable("ContaminantInMylar", srim_p1c_mylar.Data());
    } else {
        std::cerr << srim_p1c_mylar <<" not found!\n";
        return;
    }
    TString srim_p1c_IC_gas = TString::Format("../Calibrations/SRIM/17O_CF4_%0.fmbar.txt", pressureIC);
    if (std::ifstream(srim_p1c_IC_gas.Data())) {
        srim->ReadTable("ContaminantInCF4", srim_p1c_IC_gas.Data());
    } else {
        std::cerr << srim_p1c_IC_gas <<" not found!\n";
        return;
    }    


    //---- SIMULATION STARTS HERE
    ROOT::EnableImplicitMT();

    // timer
    TStopwatch timer {};
    timer.Start();
    // print fancy info
    std::cout << BOLDMAGENTA << "Running for Ex = " << Ex << " MeV" << RESET << '\n';
    std::cout << BOLDGREEN;
    const int percentPrint {5};
    int step {iterations / (100 / percentPrint)};
    int nextPrint {step};
    int percent {};


    //Histograms
        auto hT1Initial {HistConfig::dE.GetHistogram()};
        auto hTCInitial {HistConfig::dE.GetHistogram()};
        auto hT1AfterIC {HistConfig::dE.GetHistogram()};
        auto hTCAfterIC {HistConfig::dE.GetHistogram()};

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

        // Beam according to its sigma
        auto T1Initial {runner.RandomizeBeamEnergy(
                T1 * p1.GetAMU(),
                sigmaPercentBeam * T1 * p1.GetAMU())}; // T1 in Mev / u * mass of beam in u = total kinetic energy
        // std::cout << "T1Initial: " << T1Initial/p1.GetAMU() << '\n';
                
        hT1Initial->Fill(T1Initial/p1.GetAMU());
        
        // Pass beam through IC
        auto T1AfterIC {dEthroughIC(p1, T1Initial, stragglingInIC, srim)};

        hT1AfterIC->Fill(T1AfterIC/p1.GetAMU());

        // If we want to track the contaminant 17O
        if (trackContamination){
            auto TCInitial {runner.RandomizeBeamEnergy(
                T1 * p1c.GetAMU(),
                sigmaPercentBeam * T1 * p1c.GetAMU())}; 

            hTCInitial->Fill(TCInitial/p1c.GetAMU());
        
            auto TCAfterIC {dEthroughIC(p1c, TCInitial, stragglingInIC, srim)};
            hTCAfterIC->Fill(TCAfterIC/p1c.GetAMU());

        }

        // std::cout << "T1AfterIC: " << T1AfterIC/p1.GetAMU() << '\n';

    

    }
  

    //plotting
    auto* c0 {new TCanvas("c0", "Canvas for Energy Loss Plots")};

    hT1Initial->SetLineColor(kBlue);
    hT1Initial->SetStats(0);
    auto* temp1 = hT1Initial->DrawClone(); // keep clone in temp to reuse in the legend
    hT1AfterIC->SetLineColor(kRed);
    auto* temp2 = hT1AfterIC->DrawClone("same");

    auto* l0 = new TLegend(0.1,0.7,0.5,0.9); 
    l0->AddEntry(temp1, "Initial Beam Energy","l");
    l0->AddEntry(temp2, "EBeam after IC","l");

    if (trackContamination){
        hTCInitial->SetLineColor(kOrange);
        auto* temp3 = hTCInitial->DrawClone("same");
        hTCAfterIC->SetLineColor(kMagenta);
        auto* temp4 = hTCAfterIC->DrawClone("same");
        l0->AddEntry(temp3, "17O Initial Energy","l");
        l0->AddEntry(temp4, "17O after IC","l");
    }

    l0->Draw();


    timer.Stop();
    timer.Print();

}