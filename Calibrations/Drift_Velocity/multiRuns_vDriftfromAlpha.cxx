#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

#include "./vDriftfromAlpha.cxx"

void multiRuns_vDriftfromAlpha()
{
    /* List of runs for drift measurements
    run 	Vdt (V) 	Vdb (V) 	Vm (V)                  comment
    17 	    5500 	    430 	    430     700 mbar
    18 	    5800 	    430 	    430     700 mbar
    19 	    6100 	    430 	    430     700 mbar        1 CoBo failed during the run
    20 	    6100 	    430 	    430     700 mbar
    21 	    5950 	    430 	    430     700 mbar
    62      6450        480         480     755 mbar
    */

    std::map<int, double> data = {{17, (5500 - 430)}, {18, (5800 - 430)}, {20, (6100 - 430)}, {21, (5950 - 430)}};

    std::vector<double> drifts;
    auto* hDrifts = new TGraph();
    auto nPoints {0};

    for(const auto& m : data)
    {
        auto drift = vDriftfromAlpha(m.first, false);
        drifts.push_back(drift);
        std::cout << m.second << " " << drift << std::endl;
        hDrifts->SetPoint(nPoints, m.second, drift.first);
        hDrifts->SetPointError(nPoints, 0, drift.second);
        nPoints++;
    }

    auto c3 = new TCanvas("c3", "Drifts vs E", 1400, 800);
    c3->cd();

    hDrifts->SetMarkerStyle(20);
    hDrifts->SetMarkerSize(1.2);
    hDrifts->SetMarkerColor(kBlue);
    hDrifts->SetTitle("Drift vs E;E [V/cm];Drift [cm/us]");

    hDrifts->Draw("AP");
    c3->Update();
    c3->SaveAs("vdrift_runs17-21.png");

     // Compare with graph from Garfield
    std::ifstream fIn("./Inputs/H2-95_iC4H10-5_950mbar.print");
    if(!fIn.is_open())
    {
        std::cout << "Did not find Garfield file" << std::endl;
        return;
    }
    auto* hTheo = new TGraph();
    int n {0};
    std::string line;
    int iline {0};
    while(std::getline(fIn, line))
    {
        iline++;
        if(iline < 13)
            continue;
        if(iline > 32)
            break;
        std::istringstream iss(line);
        double E, V;
        if(!(iss >> E >> V))
            continue;
        hTheo->SetPoint(n, E, V);
        n++;
    }
    fIn.close();
    c3->cd();
    hTheo->SetLineColor(kRed);
    hTheo->Draw("L same");

}