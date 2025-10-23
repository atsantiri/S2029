// comparing calibrations from run 1 and 66

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

std::pair<std::vector<double>, std::vector<double>> readFile(TString fIn)
{
    std::ifstream file(fIn);
    if(!file.is_open())
    {
        std::cerr << "Could not open " << fIn << "\n";
        return {};
    }

    std::vector<double> inter, sl;
    std::string name;
    double n1, n2;
    int line = 0;

    while(file >> name >> n1 >> n2)
    {
        if(line % 2 == 0)
        { // _E lines
            inter.push_back(n1);
            sl.push_back(n2);
        }
        ++line;
    }

    return {inter, sl};
}

int comparison(const char* det_type)
{

    auto [int1, sl1] = readFile(TString::Format("Outputs/s2029_%s_redo.dat", det_type));
    auto [int2, sl2] = readFile(TString::Format("Outputs/s2029_%s_run66.dat", det_type));

    if((int1.empty() && sl1.empty()) || (int2.empty() && sl2.empty()))
        return -1;

    std::vector<double> slDiff, intDiff;
    for(size_t i = 0; i < int1.size(); ++i)
    {
        // std::cout << i << ": i1= " << int1[i] << ", s1= " << sl1[i] << "\n";


        double inter {(int1[i] - int2[i]) * 200 / (int1[i] + int2[i])};
        if(abs(inter) <= 100)
            intDiff.push_back(inter);
        else{
            intDiff.push_back(0);
            std::cout<<"Det "<<det_type<<"_"<<i<<" large intercept diff"<<std::endl;
        }


        double slope {(sl1[i] - sl2[i]) * 200 / (sl1[i] + sl2[i])};
        if(abs(slope) <= 100)
            slDiff.push_back(slope);
        else{
            slDiff.push_back(0);
            std::cout<<"Det "<<det_type<<"_"<<i<<" large slope diff"<<std::endl;
        }
    }

    TCanvas* c = new TCanvas(TString::Format("c%s", det_type), det_type, 800, 400);
    c->DivideSquare(2);
    std::vector<double> det(intDiff.size());
    for(size_t i = 0; i < intDiff.size(); ++i)
        det[i] = i;


    c->cd(1);
    TGraph* gSl = new TGraph(det.size(), det.data(), slDiff.data());
    gSl->SetTitle("Relative Slope Diff;Detector # ;Diff %");
    gSl->SetMarkerStyle(20);
    gSl->SetLineColor(kBlue);
    gSl->GetYaxis()->SetRangeUser(*std::min_element(slDiff.begin(), slDiff.end()) * 0.8,
                                  *std::max_element(slDiff.begin(), slDiff.end()) * 1.1);
    gSl->Draw("APL");

    c->cd(2);
    TGraph* gInt = new TGraph(det.size(), det.data(), intDiff.data());
    gInt->SetTitle("Relative Intercept Diff;Detector # ;Diff %");
    gInt->SetMarkerStyle(21);
    gInt->SetLineColor(kRed);
    gInt->GetYaxis()->SetRangeUser(*std::min_element(intDiff.begin(), intDiff.end()) * .8,
                                   *std::max_element(intDiff.begin(), intDiff.end()) * 1.1);
    gInt->Draw("APL");


    return 0;
}

int compareCalibrations()
{
    int ret = 0;
    if(comparison("f0") != 0)
        ret |= 1 << 0; // bit 0 = f0
    if(comparison("l0") != 0)
        ret |= 1 << 1; // bit 1 = l0
    if(comparison("r0") != 0)
        ret |= 1 << 2; // bit 2 = r0
    return ret;
}
