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
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

#include <fstream>

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

    std::vector<std::array<double, 3>> data = {{17, (5500. - 430.) / 23.5, 700},
                                               {18, (5800. - 430.) / 23.5, 700},
                                               {20, (6100. - 430.) / 23.5, 700},
                                               {21, (5950. - 430.) / 23.5, 700},
                                               {62, (6450. - 480.) / 23.5, 760}};

    double normP = 760.; // runs are with different pressures than garfield so everything has to be scaled

    std::vector<std::pair<double, double>> drifts;
    auto* hDrifts = new TGraphErrors();
    auto nPoints {0};

    for(const auto& m : data)
    {
        auto drift = vDriftfromAlpha(m[0], false);
        drifts.push_back(drift);
        std::cout << "For field E: " << m[1] << ", drift is: " << drift.first << " +- " << drift.second << " cm/us"
                  << std::endl;
        hDrifts->SetPoint(nPoints, m[1] * normP / m[2], drift.first);
        hDrifts->SetPointError(nPoints, 0, drift.second);
        nPoints++;
    }

    auto c3 = new TCanvas("c3", "Drifts vs E", 1400, 800);
    c3->cd();

    hDrifts->SetMarkerStyle(20);
    hDrifts->SetMarkerSize(1.2);
    hDrifts->SetMarkerColor(kBlue);

    // Garfield file
    std::ifstream fIn("./Inputs/Garfield_H2-iC4H10_95-5_618mbar.txt");
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
        hTheo->SetPoint(n, E * normP / 618., V);
        // std::cout << E << " " << V << std::endl;
        n++;
    }
    fIn.close();
    c3->cd();
    hTheo->SetLineColor(kRed);
    hTheo->SetTitle("Drift vs E;E [V/cm];Drift [cm/us]");
    hTheo->GetXaxis()->SetRangeUser(100,400);
    hTheo->GetYaxis()->SetRangeUser(.4,1.2);
    hTheo->Draw("AL");

    hDrifts->Draw("P same");
    c3->Update();
    c3->SaveAs("vdrift_all.png");
}