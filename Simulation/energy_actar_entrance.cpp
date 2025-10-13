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

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStopwatch.h>

#include <ActRunner.h>

#include <fstream>

#include "../PostAnalysis/HistConfig.h"


void energy_actar_entrance()
{

    // Set number of iterations
    const int iterations {static_cast<int>(1e6)};
    gRandom->SetSeed();
    // Runner: contains utility functions to execute multiple actions
    ActSim::Runner runner(nullptr, nullptr, gRandom, 0);

    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};

    ActPhysics::Particle p {"17F"};
    ActPhysics::Particle p2 {"1H"};

    // // Load SRIM Tables
    auto* srim {new ActPhysics::SRIM()};

    srim->ReadTable("mixture", "SRIM/17F_H2-iC4H10_95-5_740mbar.txt");


    //---- SIMULATION STARTS HERE
    // ROOT::EnableImplicitMT(16); // big for loops dont run in parallel

    // timer
    TStopwatch timer {};
    timer.Start();
    // print fancy info
    std::cout << BOLDGREEN;
    const int percentPrint {5};
    int step {iterations / (100 / percentPrint)};
    int nextPrint {step};
    int percent {};

    auto hEcm_dist {HistConfig::Ecm_dist.GetHistogram()};
    auto xmin {256};
    auto xmax {0};
    auto TInitial {srim->EvalEnergy("mixture", 256)/p.GetAMU()};


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
        auto Trandom {runner.RandomizeBeamEnergy(TInitial * p.GetAMU(), .005 * TInitial * p.GetAMU())}; // T1 in Mev / u * mass of beam in u = total kinetic energy


        // Find the range in the active are where the 6.15 MeV resonance is reached, should be Ecm = 2.23 MeV -> Tbeam
        // = 2.36 MeV
        int nsteps = 250;
        for(int i = 0; i < nsteps; i++)
        {
            double x = (i + 1) / static_cast<double>(nsteps) * 256; // distance travelled in active area in mm
            auto currentT {srim->Slow("mixture", Trandom, x)};
            auto Ecm = currentT / p.GetAMU() * (p.GetAMU() * p2.GetAMU()) / (p.GetAMU() + p2.GetAMU());
            hEcm_dist->Fill(x, Ecm);
            if(Ecm < 2.24 && Ecm > 2.22)
            {
                if(x < xmin)
                {
                    xmin = x;
                }
                if(x > xmax)
                {
                    xmax = x;
                }
            }
        }
    }
    //    fOut.close();

    // plotting
    std::cout << "\n Resonance region: " << xmin << " to " << xmax << " mm" << std::endl;
    auto* c1 {new TCanvas("c1", "Energy in active area")};
    hEcm_dist->SetStats(0);
    hEcm_dist->DrawClone("COLZ");
    TLine* resonance = new TLine(hEcm_dist->GetXaxis()->GetXmin(), 2.23, hEcm_dist->GetXaxis()->GetXmax(), 2.23);
    resonance->SetLineColor(kRed);
    resonance->SetLineWidth(2);
    resonance->SetLineStyle(2);
    resonance->Draw("same");
    c1->SetLogz();
    timer.Stop();
    timer.Print();
}