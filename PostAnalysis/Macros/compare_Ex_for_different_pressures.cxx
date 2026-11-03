#include "ActMergerData.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TH2.h"
#include "TLegend.h"

#include "../HistConfig.h"

std::pair<double, double> getSigma(TH1* hist)
{
    std::pair<double, double> ret;
    TF1* fE = new TF1("fE", "gaus(0)", -0.1, .2);
    hist->Fit(fE, "QRN");
    ret.first = fE->GetParameter(1);
    ret.second = fE->GetParameter(2);
    return ret;
}

void compare_Ex_for_different_pressures()
{
    // For 770 mbar-> EbeamIni = 3.89+-0.03 MeV/u
    // For 765 mbar-> EbeamIni = 3.87+-0.03 MeV/u
    // For 760 mbar-> EbeamIni = 3.85+-0.03 MeV
    // For 755 mbar-> EbeamIni = 3.83+-0.03 MeV
    // For 740 mbar-> EbeamIni = 3.78+-0.03 MeV
    // For 730 mbar-> EbeamIni = 3.75+-0.03 MeV
    std::vector<std::pair<double, double>> set {
        {3.75, 730}, {3.78, 740}, {3.83, 755}, {3.85, 760}, 
    {3.89, 770}};

    std::vector<int> colors = {616, 616, 800 - 3, 632, 432+1};

    ROOT::EnableImplicitMT();

    auto* c {new TCanvas("c", "", 1000, 800)};
    auto leg = new TLegend(0.6, 0.7, 0.88, 0.88);

    for(int i = 0; i < set.size(); i++)
    {
        ROOT::RDataFrame df {"Final_Tree", TString::Format("../Outputs/tree_ex_17F_p_p_%.2f_sil.root", set[i].first)};
        auto hEx {df.Filter([](ActRoot::MergerData& m) { return m.fLight.IsFilled() == true; }, {"MergerData"})
                      .Histo1D(HistConfig::Ex, "RecEx")};
        auto hc {(TH1*)hEx->Clone()};

        if(i == 0)
        {
            hc->Draw();
            leg->AddEntry(hc, TString::Format("%.0f mbar", set[i].second), "l");
        }
        else
        {
            hc->SetLineColor(colors[i]);
            hc->Draw("same");
            leg->AddEntry(hc, TString::Format("%.0f mbar", set[i].second), "l");
        }
        auto sig {getSigma(hc)};
        std::cout << "P = " << set[i].second << "mbar: FWHM = " << 2.355 * sig.second << std::endl;
    }
    leg->Draw();
}