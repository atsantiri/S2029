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

std::pair<std::vector<double>, std::vector<double>> readRes(TString fIn)
{
    std::ifstream file(fIn);
    if(!file.is_open())
    {
        std::cerr << "Could not open " << fIn << "\n";
        return {};
    }

    std::vector<double> res, err;
    std::string name, temp1, temp2;
    double n1, n2;
    int line = 0;

    while(file >> name >> n1 >> temp1 >> n2 >> temp2)
    {
        res.push_back(n1);
        err.push_back(n2);
        ++line;
    }

    return {res, err};
}


int comparison(const char* det_type)
{

    auto [int1, sl1] = readFile(TString::Format("Outputs/s2029_%s_run1.dat", det_type));
    auto [int2, sl2] = readFile(TString::Format("Outputs/s2029_%s_run66.dat", det_type));

    auto [res1, err1] = readRes(TString::Format("Outputs/resolution_%s_run1.dat", det_type));
    auto [res2, err2] = readRes(TString::Format("Outputs/resolution_%s_run66.dat", det_type));


    if((int1.empty() && sl1.empty()) || (int2.empty() && sl2.empty()) || (res1.empty() && err1.empty()) ||
       (res2.empty() && err2.empty()))
        return -1;

    std::vector<double> slDiff, intDiff;
    for(size_t i = 0; i < int1.size(); ++i)
    {
        // std::cout << i << ": i1= " << int1[i] << ", s1= " << sl1[i] << "\n";
        // std::cout << i << ": i1= " << res1[i] << ", s1= " << err1[i] << "\n";


        double inter {(int1[i] - int2[i]) * 200 / (int1[i] + int2[i])};
        if(abs(inter) <= 100)
            intDiff.push_back(inter);
        else
        {
            intDiff.push_back(0);
            // std::cout<<"Det "<<det_type<<"_"<<i<<" large intercept diff "<<inter<<"%"<<std::endl;
        }


        double slope {(sl1[i] - sl2[i]) * 200 / (sl1[i] + sl2[i])};
        if(abs(slope) <= 100)
            slDiff.push_back(slope);
        else
        {
            slDiff.push_back(0);
            std::cout << "Det " << det_type << "_" << i << " large slope diff " << slope << "%" << std::endl;
        }

        double res_diff {abs(res1[i] - res2[i])};
        if(abs(res_diff) >= 5)
        {
            std::cout << "Det " << det_type << "_" << i << " large resolution diff "<< res1[i]<<" vs "<<res2[i] << std::endl;
        }

        if(res1[i] > 0 && res2[i] > 0 && res1[i] < res2[i])
            std::cout << "For Det " << det_type << "_" << i << " use run 1" << std::endl;
        else if(res1[i] > 0 && res2[i] > 0 && res1[i] > res2[i])
            std::cout << "For Det " << det_type << "_" << i << " use run 66" << std::endl;
        else
            std::cout << "For Det " << det_type << "_" << i << " I'm confused" << std::endl;

    }


    return 0;
}

int compareCalibrations()
{
    int ret = 0;
    std::cout<<"================ F0 ================"<<std::endl;
    if(comparison("f0") != 0)
        ret |= 1 << 0; // bit 0 = f0
    std::cout<<"================ L0 ================"<<std::endl;
    if(comparison("l0") != 0)
        ret |= 1 << 1; // bit 1 = l0
    std::cout<<"================ R0 ================"<<std::endl;
    if(comparison("r0") != 0)
        ret |= 1 << 2; // bit 2 = r0
    return ret;
}
