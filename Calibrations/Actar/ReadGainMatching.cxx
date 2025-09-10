#include "ActCalibrationManager.h"
#include "ActTPCLegacyData.h"
#include "ActTPCParameters.h"

#include "TCanvas.h"
#include "TChain.h"
#include "TH2.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>
#include <iomanip>

void FillHistogram(ActRoot::CalibrationManager* calman, ActRoot::TPCParameters* pars, MEventReduced* evt, TH2D* h,
                   bool isMatched)
{

    // iterate over hits
    for(int it = 0, size = evt->CoboAsad.size(); it < size; it++)
    {
        // Global Channel ID: [ CoBo (5 bits: 0-15) | AsAd (2 bits: 0-3) | AGET (2 bits: 0-3) | Channel (7 bits: 0-127, but only 64+4 are real) ]
        int co = evt->CoboAsad[it].globalchannelid >> 11;                                   // Shift right 11 bits -> highest 5 bits: CoBo
        int as = (evt->CoboAsad[it].globalchannelid - (co << 11)) >> 9;                     // Subtract CoBo and shift right 9 bits. The next 2 bits is AsAd.
        int ag = (evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9)) >> 7;         // Subtract CoBo and AsAd, then shift right 7 bits. Next 2 bits is AGET.
        int ch = evt->CoboAsad[it].globalchannelid - (co << 11) - (as << 9) - (ag << 7);    // Remove all of the above and the channel number remains.
        int where = co * pars->GetNBASAD() * pars->GetNBAGET() * pars->GetNBCHANNEL() +
                    as * pars->GetNBAGET() * pars->GetNBCHANNEL() + ag * pars->GetNBCHANNEL() + ch;

        if((co != 31) && (co != 16))
        {
            auto xval {calman->ApplyLookUp(where, 4)};
            auto yval {calman->ApplyLookUp(where, 5)};
            for(int hit = 0, otherSize = evt->CoboAsad[it].peakheight.size(); hit < otherSize; hit++)
            {
                if((yval != -1) && (xval != -1))
                {
                    double z_position {evt->CoboAsad[it].peaktime[hit]};
                    if(z_position > 0.)
                    {
                        auto Qiaux {evt->CoboAsad[it].peakheight[hit]};
                        // Fill histogram
                        if(isMatched)
                            Qiaux = calman->ApplyPadAlignment(where, Qiaux);
                        // Fill histogram

                        if(Qiaux >= 500) // delete baseline
                        {
                            h->Fill(where, Qiaux);
                        }
                    }
                }
            }
        }
    }
}

void ReadGainMatching(bool isMatched = true)
{

    // Get the data into TChain
    auto chain {new TChain("ACTAR_TTree")};
    std::vector<int> runs {12, 13, 14, 15};
    for(const auto& run : runs)
    {
        chain->Add(TString::Format("../../RootFiles/Raw/Tree_Run_%04d_Merged.root", run));
    }

    // Set the parameters
    ActRoot::TPCParameters tpc {"Actar"};
    // Calibration manager
    ActRoot::CalibrationManager calman {};
    calman.ReadLookUpTable("./LT.txt"); // See README for format of LT table
    if(isMatched)
        calman.ReadPadAlign("./Outputs/gain_matching_s2029.txt");

    // Create histogram
    auto* h {new TH2D {"h", "pads;Channel number;", 17408, 0, 17408, 1000, 0, 5000}};
    if (isMatched)
        h->SetYTitle("Calibrated charge [u.a.]");
    else
        h->SetYTitle("Uncalibrated charge [u.a.]");

    // Set MEventReduced
    MEventReduced* evt {new MEventReduced};
    chain->SetBranchAddress("data", &evt); // get the column from the chain that we gonna get data from

    for(long int entry = 0, maxEntry = chain->GetEntries(); entry < maxEntry; entry++)
    {
        std::cout << "\r"
                  << "At percent : " << std::fixed << std::setprecision(1)<< ((double)entry / maxEntry) * 100 << " %" << std::flush;
        chain->GetEntry(entry); // get  the data from the chain and write it to evt variable
        FillHistogram(&calman, &tpc, evt, h, isMatched);
    }

    // Plot
    gStyle->SetOptStat(0);
    auto* c0 {new TCanvas {"c0", "Gain matching canvas"}};
    h->Draw("colz");

    if(!isMatched)
        h->SaveAs("./Inputs/gain.root");
    else
        h->SaveAs("./Outputs/aligned.root");
}
