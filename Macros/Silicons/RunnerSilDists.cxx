#include "ActColors.h"

#include "TROOT.h"
#include "TString.h"

#include <iostream>
#include <string>


void RunnerSilDists(TString what = "")
{
    std::string detector {"r0"};

    std::cout << BOLDGREEN << "···· RunnerSilDists ····" << '\n';
    std::cout << "-> Silicon  : " << detector << '\n';
    std::cout << "-> What     : " << what << '\n';
    std::cout << "······························" << RESET << '\n';

    auto args {TString::Format("(\"%s\")", detector.c_str())};
    TString func {};
    TString ext {".cxx"};

    // DistRun: Creates histograms for different sil centres
    if(what.Contains("1"))
    {
        func = "DistRun";
        gROOT->LoadMacro(func + ext);
        gROOT->ProcessLine(func + args);
    }
    // DistPlot: Fits the x/y and z projections 
    if(what.Contains("2"))
    {
        func = "DistPlot";
        gROOT->LoadMacro(func + ext);
        gROOT->ProcessLine(func + args);
    }
    // DistDebug: From the heights/widths of the fits it finds which one is closest to the real value
    if(what.Contains("3"))
    {
        func = "DistDebug";
        gROOT->LoadMacro(func + ext);
        gROOT->ProcessLine(func + args);
    }
}
