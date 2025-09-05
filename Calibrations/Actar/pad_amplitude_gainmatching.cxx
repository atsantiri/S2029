#include "ActDataManager.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include "ROOT/RDataFrame.hxx"
// #include "ROOT/TThreadedObject.hxx"

#include <vector>
#include "TH1D.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TF1.h"
#include <map>
#include <fstream>


void pad_amplitude_gainmatching()
{
    auto* c0 {new TCanvas {"c0", "Pad amplitude gainmatching"}};
    c0->DivideSquare(5);

    int p {1};
    // make a map of TGraphs, one per pad, to keep the amplitude vs reference amplitude for the fits
    std::map<int, TGraph*> padGraphs;
    for (int pad = 0; pad < 16384; ++pad) {
        padGraphs[pad] = new TGraph();
    }
    int referencePad = 9453; 
    int testPad {5};

    for (int run=12; run<=16; ++run) {
        ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EReadTPC};
        dataman.SetRuns(run,run);
        auto chain {dataman.GetChain()};
        ROOT::RDataFrame df {*chain};

        // auto* hq {new TH1D {"hq", "pulser run;pad ID ; Qave", 16384, 0, 16383}};
        auto hq {new TProfile(TString::Format("hq%d", run), TString::Format("Pulser run %d;pad ID;Average Q(ch)",run), 16384, 0, 16383)}; //TProfile computes avrages instead of sum

        auto df2 {df.Define("padIDs", 
                    [](const ActRoot::TPCData& tpc) {
                        std::vector<int> pads;
                        for (const auto& v : tpc.fRaw)
                            pads.push_back((int)v.GetPosition().X() * 128 + (int)v.GetPosition().Y());
                        return pads;
                    }, 
                    {"TPCData"})
                .Define("charges", 
                    [](const ActRoot::TPCData& tpc) {
                        std::vector<double> q;
                        for (const auto& v : tpc.fRaw)
                            q.push_back(v.GetCharge());
                        return q;
                    }, 
                    {"TPCData"})};

        df2.Foreach(
            [&](const std::vector<int>& pads, const std::vector<double>& charges){
                for (size_t i=0; i<pads.size(); ++i)
                    hq->Fill(pads[i], charges[i]);  // TProfile stores average per bin
            },
            {"padIDs","charges"}
        );

        c0->cd(p++);
        hq->Draw();
        hq->GetYaxis()->SetRangeUser(0, 200);

        double refCharge = hq->GetBinContent(referencePad+1);
        for (int pad = 0; pad < 16384; ++pad) {
            if (hq->GetBinEntries(pad+1) > 0) { 
                double thisCharge = hq->GetBinContent(pad+1);
                padGraphs[pad]->AddPoint(refCharge, thisCharge);
            }
        }

        c0->cd(6);
        padGraphs[testPad]->SetMarkerStyle(5);
        padGraphs[testPad]->SetMarkerSize(3);
        padGraphs[testPad]->SetTitle(TString::Format("Test pad %d amplitude vs ref pad %d", testPad, referencePad));
        padGraphs[testPad]->GetXaxis()->SetTitle(TString::Format("Ref pad %d ", referencePad));
        padGraphs[testPad]->GetYaxis()->SetTitle(TString::Format("Test pad %d ", testPad));
        if (run == 12) 
            padGraphs[testPad]->Draw("AP");
        else
            padGraphs[testPad]->Draw("P same");

        // Find which pad to normalize to based on the overall average charge
        double sum = 0;
        int n = 0;
        for (int bin = 1; bin <= hq->GetNbinsX(); ++bin) { 
            double val = hq->GetBinContent(bin);
            if (hq->GetBinEntries(bin) > 0) { 
                sum += val;
                ++n;
            }
        }

        double overallQave = (n > 0) ? sum / n : 0;
        int refPad = -1;
        double minDiff = 1e30;
        for (int b = 1; b <= hq->GetNbinsX(); ++b) {
            if (hq->GetBinEntries(b) > 0) {
                double val = hq->GetBinContent(b);
                double diff = std::abs(val - overallQave);
                if (diff < minDiff) {
                    minDiff = diff;
                    refPad = b - 1; 
                }
            }
        }
        std::cout<< "Run " << run << ": "
                << "Overall average: " << overallQave
                << ", closest pad ID: " << refPad
                << ", pad value: " << hq->GetBinContent(refPad+1) 
                << ", pad 9453: " << hq->GetBinContent(9453+1) << std::endl;

        // Will use pad 9453 as reference for now (it is close to the average for all runs)
    }

    // Test pad fit
    c0->cd(6);
    TF1 *f = new TF1("lin","[0]+[1]*x");
    padGraphs[testPad]->Fit(f, "Q"); // quiet
    double offset = f->GetParameter(0);
    double slope  = f->GetParameter(1);
    std::cout << "Test pad " << testPad << 
                ": offset=" << offset
                << ", slope=" << slope << std::endl;


    // Fit all pads
    std::ofstream streamer {"./pad_alignment_s2029.txt"};
    for (auto& [pad, gr] : padGraphs) {
        if (gr->GetN() > 3) { 
            TF1 *f = new TF1("lin","[0]+[1]*x");
            gr->Fit(f, "Q"); // quiet
            double offset = f->GetParameter(0);
            double slope  = f->GetParameter(1);
            streamer << slope << "\t" << offset << std::endl;
        }
        else 
            streamer << "0 \t 0" << std::endl;
    }
    streamer.close();

}
