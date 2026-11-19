#include "ActMergerData.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/TThreadedObject.hxx>

#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"

#include "../HistConfig.h"
#include "../Settings.h"

std::pair<double, double> getSigma(TH1* hist)
{
    std::pair<double, double> ret;
    TF1* fE = new TF1("fE", "gaus(0)", -0.7, .7);
    hist->Fit(fE, "QRN");
    ret.first = fE->GetParameter(1);
    ret.second = fE->GetParameter(2);
    return ret;
}

void debugF0dist()
{
    std::vector<int> colors = {TColor::GetColor("#08306B"), TColor::GetColor("#08519C"), TColor::GetColor("#2171B5"),
                               TColor::GetColor("#4292C6"), TColor::GetColor("#6BAED6"), TColor::GetColor("#9ECAE1"),
                               TColor::GetColor("#C6DBEF")};
    std::vector<int> set {300, 310, 320, 330, 340, 360, 370};
    ROOT::EnableImplicitMT();

    auto* c0 {new TCanvas("c0", "Ex", 1000, 800)};
    auto* c1 {new TCanvas("c1", "2D plots", 1000, 800)};
    auto leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    c1->DivideSquare(set.size());
    std::cout << "F0 pos \t FWHM \tcentroid " << std::endl;
    for(int i = 0; i < set.size(); i++)
    {
        ROOT::RDataFrame df {"Final_Tree", TString::Format("../Outputs/tree_ex_17F_p_p_3.87_sil_%dF0.root", set[i])};
        auto hEx {df.Filter([](ActRoot::MergerData& m)
                            { return (m.fLight.GetNLayers() == 1 && m.fLight.GetLayer(0) == "f0"); }, {"MergerData"})
                      .Histo1D(HistConfig::Ex, "RecEx")};
        c0->cd();
        auto hc {(TH1*)hEx->Clone()};

        if(i == 0)
        {
            hc->GetYaxis()->SetRangeUser(0, hc->GetMaximum() * 1.1);
            hc->Draw();
            leg->AddEntry(hc, TString::Format("%d", set[i]), "l");
        }
        else
        {
            hc->SetLineColor(colors[i]);
            hc->Draw("same");
            leg->AddEntry(hc, TString::Format("%d", set[i]), "l");
        }

        auto sig {getSigma(hc)};
        std::cout << set[i] << std::setprecision(3) << "\t" << 2.355 * sig.second << "\t" << sig.first << std::endl;


        auto h2D {df.Filter([](ActRoot::MergerData& m) { return (m.fLight.GetNLayers() == 1); }, {"MergerData"})
                      .Histo2D(HistConfig::ExRPx, "fRP.fCoordinates.fX", "RecEx")};
        h2D->SetTitle(TString::Format("F0 @ %d", set[i]));
        h2D->GetYaxis()->SetRangeUser(-3., 3.);
        c1->cd(i + 1);
        h2D->DrawClone();
    }
    c0->cd();
    leg->Draw();
}