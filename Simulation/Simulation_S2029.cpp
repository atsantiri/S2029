// Based Simulation_E796.cpp

#include "ActColors.h"
#include "ActCrossSection.h"
#include "ActDecayGenerator.h"
#include "ActKinematicGenerator.h"
#include "ActKinematics.h"
#include "ActLine.h"
#include "ActParticle.h"
#include "ActSRIM.h"
#include "ActSilData.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActUtils.h"
#include <ActRunner.h>

#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TFile.h>
#include <fstream>
#include <TStopwatch.h>
#include <TH2D.h>
#include <TLine.h>

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
    auto eneAfterMylar {srim->Slow(tableKey1.Data(), energy, 11.5e-3)}; // 11.5 um mylar in IC 
    // 36 mm CF4 gas
    auto eneAfterIC {srim->Slow(tableKey2.Data(), eneAfterMylar, 36)};
    
    return eneAfterIC;
}

double dEthroughCFA(const ActPhysics::Particle& particle, double energy, bool applyStraggling, ActPhysics::SRIM* srim) {
    if (!applyStraggling) return energy;

    // Choose the right SRIM table key based on particle name
    TString tableKey1, tableKey2;
    if (particle.GetName() == "17F" || particle.GetName() == "20Mg") {
        tableKey1 = "beamInMylar";
        tableKey2 = "beamIniC4H10";
    } else if (particle.GetName() == "17O" || particle.GetName() == "20Na") {
        tableKey1 = "ContaminantInMylar";
        tableKey2 = "ContaminantIniC4H10";
    } else {
        std::cerr << "Unknown particle in dEthroughIC: " << particle.GetName() << '\n';
        return energy; 
    }
    
    // 4.8 um mylar
    auto eneAfterMylar {srim->Slow(tableKey1.Data(), energy, 4.8e-3)}; 
    // 9.6 mm iC4H10 (isobutane) gas
    auto eneAfterCFA {srim->Slow(tableKey2.Data(), eneAfterMylar, 9.6)};
    
    return eneAfterCFA;
}



void Simulation_S2029(const std::string& beam = "17F", double T1 = 5.5, double Ex = 0){

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
    double pressureACTAR {600}; // Pressure in mbar    
    bool trackContamination {false}; // If true the 17O contaminant will be propagated along with the 17F primary beam

    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    // auto* sils {new ActPhysics::SilSpecs};
    // sils->ReadFile("../configs/silspecs.conf");
    // sils->Print();

    // Resolutions
    const double sigmaSil {0.060 / 2.355};
    const double sigmaPercentBeam {0.005};
    const double sigmaAngleLight {0.95 / 2.355};

    // Kinematics
    ActPhysics::Particle p1 {beam};
    ActPhysics::Particle p1c {"17O"}; // Contaminant 
    ActPhysics::Particle p2 {"1H"};
    // ActPhysics::Particle p3 {arglight};
    // // Automatically compute 4th particle
    // ActPhysics::Kinematics kaux {p1, p2, p3};
    // ActPhysics::Particle p4 {kaux.GetParticle(4)};
    // // Binary kinematics generator
    // ActSim::KinematicGenerator kingen {
    //     p1, p2, p3, p4, (protonPS > 0 ? protonPS : 0), (neutronPS > 0 ? neutronPS : (pdphase || dtphase ? 1 : 0))};
    // kingen.Print();


    // // Load SRIM Tables
    auto* srim {new ActPhysics::SRIM()};

    auto load_table = [&](const char* label, const TString& path) {
        if (std::ifstream(path.Data())) {
            srim->ReadTable(label, path.Data());
        } else {
            std::cerr << path << " not found!\n";
            return false;
        }
        return true;
    };

    std::vector<std::pair<const char*, TString>> tables = {
        {"beamInMylar", TString::Format("./SRIM/%s_mylar.txt", beam.c_str())},
        {"beamInCF4", TString::Format("./SRIM/%s_CF4_%.0fmbar.txt", beam.c_str(), pressureIC)},
        {"beamIniC4H10", TString::Format("./SRIM/%s_iC4H10_%.0fmbar.txt", beam.c_str(), pressureCFA)},
        {"beamInACTARgas", TString::Format("./SRIM/%s_H2-iC4H10_95-5_%.0fmbar.txt", beam.c_str(), pressureACTAR)},
        // {"ContaminantInMylar", TString::Format("./SRIM/%s_mylar.txt", p1c.GetName().c_str())},
        // {"ContaminantInCF4", TString::Format("./SRIM/%s_CF4_%.0fmbar.txt", p1c.GetName().c_str(), pressureIC)},
        // {"ContaminantIniC4H10", TString::Format("./SRIM/%s_iC4H10_%.0fmbar.txt", p1c.GetName().c_str(), pressureCFA)},
        // {"ContaminantInACTARgas", TString::Format("./SRIM/%s_H2-iC4H10_95-5_%.0fmbar.txt", p1c.GetName().c_str(), pressureACTAR)}
    };

    for (const auto& [label, path] : tables) {
        if (!load_table(label, path))
            return;
    }

    // std::ofstream fOut("dE_out.txt"); // Output test file to compare ACTPhysics Slow function with TRIM
    // fOut << "17F_IC \t 17O_IC \t 17F_CFA \t 17O_CFA \n";  


    //---- SIMULATION STARTS HERE
    //ROOT::EnableImplicitMT(16); // big for loops dont run in parallel

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
    std::map<std::string, TH1*> hSRIM;
    std::map<std::string, TH1*> hSRIMC; // to follow the contaminant
    std::vector<std::string> hSRIMLabels = {
        "TInitial",
        "TAfterIC",
        "TAfterCFA",
        "TAfterFoil",
        "TActiveAreaEntrance",
        "TActiveAreaMid",
        "TActiveAreaEnd"
    };
    auto hEcm_dist {HistConfig::Ecm_dist.GetHistogram()};
    auto xmin {256};
    auto xmax {0};

    for (const auto& label : hSRIMLabels) {
        hSRIM[label]=new TH1D(label.c_str(), label.c_str(), 200,0,5.5);
        hSRIMC[label]=new TH1D((label + "_Cont").c_str(), (label + "_Cont").c_str(), 200,0,5.5);
    }

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
                
        // hT1Initial->Fill(T1Initial/p1.GetAMU());
        hSRIM["TInitial"]->Fill(T1Initial/p1.GetAMU());

        // Pass beam through IC
        auto T1AfterIC {dEthroughIC(p1, T1Initial, stragglingInIC, srim)};
        hSRIM["TAfterIC"]->Fill(T1AfterIC/p1.GetAMU());

        // Pass beam through CFA
        auto T1AfterCFA {dEthroughCFA(p1, T1AfterIC, stragglingInCFA, srim)};
        hSRIM["TAfterCFA"]->Fill(T1AfterCFA/p1.GetAMU());

        // Entrance mylar foil of ACTAR
        auto T1AfterEntranceWindow {srim->Slow("beamInMylar", T1AfterCFA, 12e-3)}; // 12 um mylar in ACTAR entrance window
        hSRIM["TAfterFoil"]->Fill(T1AfterEntranceWindow/p1.GetAMU());

        // Space before field cage
        auto T1FieldCageEntrance {srim->Slow("beamInACTARgas", T1AfterEntranceWindow, 35)}; // 35 mm space before field cage 
        auto T1ActiveAreaEntrance {srim->Slow("beamInACTARgas", T1FieldCageEntrance, 18)}; // 18 mm space before field cage 
        hSRIM["TActiveAreaEntrance"]->Fill(T1ActiveAreaEntrance/p1.GetAMU());

        auto T1ActiveAreaMid {srim->Slow("beamInACTARgas", T1ActiveAreaEntrance, 128)};
        hSRIM["TActiveAreaMid"]->Fill(T1ActiveAreaMid/p1.GetAMU());

        auto T1ActiveAreaEnd {srim->Slow("beamInACTARgas", T1ActiveAreaEntrance, 256)}; 
        hSRIM["TActiveAreaEnd"]->Fill(T1ActiveAreaEnd/p1.GetAMU()); 

        // If we want to track the contaminant 17O
        if (trackContamination){
            auto TCInitial {runner.RandomizeBeamEnergy( T1 * p1c.GetAMU(), sigmaPercentBeam * T1 * p1c.GetAMU())}; 
            hSRIMC["TInitial"]->Fill(TCInitial/p1c.GetAMU());
        
            auto TCAfterIC {dEthroughIC(p1c, TCInitial, stragglingInIC, srim)};
            hSRIMC["TAfterIC"]->Fill(TCAfterIC/p1c.GetAMU());

            auto TCAfterCFA {dEthroughCFA(p1c, TCAfterIC, stragglingInCFA, srim)};
            hSRIMC["TAfterCFA"]->Fill(TCAfterCFA/p1c.GetAMU());

            auto TCAfterEntranceWindow {srim->Slow("ContaminantInMylar", TCAfterCFA, 12e-3)}; // 12 um mylar in ACTAR entrance window
            hSRIMC["TAfterFoil"]->Fill(TCAfterEntranceWindow/p1c.GetAMU()); 

            auto TCFieldCageEntrance {srim->Slow("ContaminantInACTARgas", TCAfterEntranceWindow, 31.6)}; // 35 mm space before field cage
            auto TCActiveAreaEntrance {srim->Slow("ContaminantInACTARgas", TCFieldCageEntrance, 20.4)}; // (12.8+7.6) mm space before field cage
            hSRIMC["TActiveAreaEntrance"]->Fill(TCActiveAreaEntrance/p1c.GetAMU());    
            auto TCActiveAreaMid {srim->Slow("ContaminantInACTARgas", TCActiveAreaEntrance, 128)};
            hSRIMC["TActiveAreaMid"]->Fill(TCActiveAreaMid/p1c.GetAMU());    
            auto TCActiveAreaEnd {srim->Slow("ContaminantInACTARgas", TCActiveAreaEntrance, 256)};
            hSRIMC["TActiveAreaEnd"]->Fill(TCActiveAreaEnd/p1c.GetAMU());
        }


        // Find the range in the active are where the 6.15 MeV resonance is reached, should be Ecm = 2.23 MeV -> Tbeam = 2.36 MeV
        int nsteps = 250;
        for (int i = 0; i < nsteps; i++) {
            double x = (i+1) / static_cast<double>(nsteps) * 256; // distance travelled in active area in mm
            auto currentT {srim->Slow("beamInACTARgas", T1ActiveAreaEntrance, x)};
            auto Ecm = currentT/p1.GetAMU()*(p1.GetAMU()*p2.GetAMU())/(p1.GetAMU() + p2.GetAMU());
            hEcm_dist->Fill(x, Ecm);
            if (Ecm<2.24 && Ecm>2.22){
                if (x<xmin){xmin=x;}
                if (x>xmax){xmax=x;}
            }
        }
    }
//    fOut.close();

    //plotting
    auto* c0 {new TCanvas("c0", "Energy Loss")};
    auto* l0 = new TLegend(0.1,0.7,0.5,0.9); 

    bool firstDraw = true;

    // Define a set of colors for the histograms
    std::vector<int> colors = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2, kBlack};
    int colorIdx = 0;

    for (const auto& label : hSRIMLabels) {
        auto* h = hSRIM[label];
        h->SetStats(0);
        int color = colors[colorIdx % colors.size()];
        h->SetLineColor(color);
        h->SetLineWidth(2);

        if (firstDraw) {
            h->SetTitle(TString::Format("%s Energy Loss", beam.c_str()));
            h->DrawClone();
            firstDraw = false;
        } else {
            h->DrawClone("same");
        }
        l0->AddEntry(h, label.c_str(), "l");

        if (trackContamination) {
            auto* hc = hSRIMC[label];
            hc->SetLineStyle(2);
            hc->SetLineColor(color);
            hc->SetLineWidth(2);
            hc->DrawClone("same");
            l0->AddEntry(hc, (label + " Cont").c_str(), "l");
        }
        ++colorIdx;
    }

    l0->Draw();
    std::cout<<"\n Resonance region: "<<xmin<<" to "<<xmax<<" mm"<<std::endl;
    auto* c1 {new TCanvas("c1", "Energy in active area")};
    hEcm_dist->SetStats(0);
    hEcm_dist->DrawClone("COLZ");
    TLine* resonance = new TLine(hEcm_dist->GetXaxis()->GetXmin(), 2.23, hEcm_dist->GetXaxis()->GetXmax(), 2.23);
    resonance->SetLineColor(kRed);
    resonance->SetLineWidth(2);
    resonance->SetLineStyle(2); 
    resonance->Draw("same");
    c1->SetLogz();
    c1->SaveAs(TString::Format("17F_energy_in_actar_%.0fmbar.png",pressureACTAR));
    timer.Stop();
    timer.Print();

}