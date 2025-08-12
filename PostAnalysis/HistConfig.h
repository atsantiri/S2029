#ifndef HistConfig_h
#define HistConfig_h

#include "ROOT/RDF/HistoModels.hxx"

#include "TString.h"

namespace HistConfig
{
using namespace ROOT::RDF;

const TH1DModel Ex {"hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (35. - (-10.)) / 300 * 1e3), 200, -10, 35};
const TH1DModel dE {"hdE", "Energy Loss;E_{x} [MeV/u];Counts", 200,0,5.5};
const TH2DModel Ecm_dist {"hEcm_dist","Energy along active area;Distance [mm];E_{CM} [MeV/u]",200,0,256,200,0,4.0};


template <typename T>
T ChangeTitle(T model, const TString& title, const TString& label = "");
} // namespace HistConfig

template <typename T>
T HistConfig::ChangeTitle(T model, const TString& title, const TString& label)
{
    auto ret {model};
    if(label.Length() > 0)
        ret.fName = model.fName + label;
    TString old {model.fTitle};
    auto end {old.First(';')};
    TString nt {title + old(end, old.Length() - end)};
    ret.fTitle = nt;
    return ret;
}

#endif 
