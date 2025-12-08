
#include "ActKinematics.h"
#include "ActMergerData.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "ROOT/RDataFrame.hxx"
#include "Rtypes.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

#include "../HistConfig.h"

double chi2 (TH1D* hRec, TH1D* hSRIM){

    double ret {0};
    // for (int i=34; i<128; i++){
    for (int i=0; i<hRec->GetNbinsX(); i++){
        double cRec = hRec->GetBinContent(i);
        double cSRIM = hSRIM->GetBinContent(i);
        if (cRec >0 && cSRIM >0)
        ret += (cRec-cSRIM)*(cRec-cSRIM)/cRec;
    }
    return ret;
}


void compare_SRIM_tables()
{

    // std::vector<double> vars {765., 760., 755., 740.};
    std::vector<double> vars {3.82};
    // std::vector<double> vars {3.70,3.71,3.72,3.73,3.74, 3.75,3.76,3.77,3.78,3.79,3.80,3.81, 3.82,3.83, 3.84, 3.85, 3.86,3.87,3.88,3.89, 3.90};
    auto* c {new TCanvas("c", "ECN", 800, 600)};

    bool first {true};
    std::vector<int> cols {2, 4, 6, 7, 8, 1, 1, 1, 1, 1, 1,1,1,1,1,1,1, 1, 1, 1};
    auto l {new TLegend(0.75, 0.65, 0.9, 0.9)};

    for(int i =0; i< vars.size(); i++)
    {
        TString fIn = TString::Format("../Outputs/tree_ex_17F_p_p_%.2f.root", vars[i]);
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame df {"Final_Tree", fIn};

        auto hEcn {df.Histo1D(HistConfig::ECN, "ECN")};
        auto hRecEcn {df.Histo1D(HistConfig::ECN, "RecECN")};
        TH1D* hrectemp = (TH1D*)hRecEcn->Clone();

        if(first)
        {
            hRecEcn->SetLineWidth(2);
            hRecEcn->DrawClone();
            first = false;
        }
        TH1D* htemp = (TH1D*)hEcn->Clone();
        std::cout << "Ebeam "<< vars[i]<<", chi2: "<<chi2(hrectemp, htemp)<<std::endl;
        htemp->SetLineColor(cols[i]);
        l->AddEntry(htemp, TString::Format("Ebeam = %.2f MeV/u", vars[i]));
        htemp->Draw("same");
    }
    l->Draw();
}