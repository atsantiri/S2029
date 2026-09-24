#include "ActDataManager.h"
#include "ActLine.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/TThreadedObject.hxx"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TPaveText.h"

#include <filesystem>

ROOT::RDataFrame getData() // If the file exists load it, otherwise create it and load it.
{
    const TString fData = "../Outputs/emittance_17F.root";
    if(std::filesystem::exists(fData.Data()))
    {
        std::cout << "Reading in file: " << fData << std::endl;
        ROOT::RDataFrame df {"Emittance_Tree", fData};
        return df;
    }
    else
    {
        std::cout << "Didn't find file: " << fData << std::endl;
        std::cout << "Processing data.conf" << std::endl;
        ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
        auto chain {dataman.GetChain()};
        auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
        chain->AddFriend(chain2.get());
        ROOT::EnableImplicitMT();

        ROOT::RDataFrame df {*chain};

        // Read conversion factors
        float padSide {2}; // mm
        ActRoot::InputParser parser {"../../configs/detector.conf"};
        auto merger {parser.GetBlock("Merger")};
        float driftFactor {static_cast<float>(merger->GetDouble("DriftFactor"))};
        std::cout << "-> DriftFactor : " << driftFactor << '\n';

        auto dff {
            df.Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 64; },
                      {"ModularData"})                       // gate on CFA_DIV
                .Filter("fClusters.fIsBeamLike.size() == 1") // gate on single track events
                .Filter("fClusters.fIsBeamLike.front() == true")
                .Define("Line",
                        [&](ActRoot::TPCData& data)
                        {
                            auto& voxels {data.fClusters.front().GetVoxels()};
                            ActRoot::Line line;
                            line.FitVoxels(voxels, true, true, true);
                            line.Scale(padSide, driftFactor);
                            return line;
                        },
                        {"TPCData"})
                .Define("Entrance", [](const ActRoot::Line& l) { return l.MoveToX(0); }, {"Line"})
                .Define("Range", "fClusters.fXRange.second") // I would in theory want to take the beamspot at the end
                                                             // of the active area, but since the beam stops in ACTAR
                                                             // want to find an average range. => Seems like the beam
                                                             // stops around 125, so I'll take the end beamspot at 100
                .Define("End", [](const ActRoot::Line& l) { return l.MoveToX(100); }, {"Line"})
                .Define("Y1", "Entrance.Y()")
                .Define("Y2", "End.Y()")
                .Define("Z1", "Entrance.Z()")
                .Define("Z2", "End.Z()")
                .Define("thetaXY",
                        [](ActRoot::Line& l)
                        {
                            auto d = l.GetDirection();
                            // atan2: XY angle, signed according to whether Z is >0 or <0>
                            return std::atan2(d.Y(), d.X()) * TMath::RadToDeg();
                        },
                        {"Line"})
                .Define("thetaXZ",
                        [](ActRoot::Line& l)
                        {
                            auto d = l.GetDirection();
                            return std::atan2(d.Z(), d.X()) * TMath::RadToDeg();
                        },
                        {"Line"})};
        ROOT::RDF::Experimental::AddProgressBar(dff);
        dff.Snapshot("Emittance_Tree", fData,
                     {"Entrance", "End", "Y1", "Y2", "Z1", "Z2", "Line", "thetaXY", "thetaXZ"});
        std::cout << "Created file: " << fData << std::endl;
        return ROOT::RDataFrame {"Emittance_Tree", fData};
    }
}


void beamEmittance()
{
    ROOT::EnableImplicitMT();
    auto df = getData();

    auto hEntry {df.Histo2D({"hEntry", "X = 0 mm;Y [mm];Z [mm]", 500, 0, 256, 300, 0, 256}, "Y1", "Z1")};
    auto hEntryZ {df.Histo1D("Z1")};
    auto hEnd {df.Histo2D({"hEnd", "X = 100 mm;Y [mm];Z [mm]", 500, 0, 256, 300, 0, 256}, "Y2", "Z2")};
    auto hThetaXZ {df.Histo1D({"hThetaXZ", "XZ Angle;;", 250, -5, 5}, "thetaXZ")};
    auto hThetaXY {df.Histo1D({"hThetaXY", "XY Angle;;", 250, -5, 5}, "thetaXY")};
    auto hYThetaXY {df.Histo2D({"hYThetaXY", "Y vs #theta_{XY};Y [mm];#theta_{XY} [#circ]", 500, 0, 256, 150, -10, 10},
                               "Y1", "thetaXY")};
    auto hYThetaXZ {df.Histo2D({"hYThetaXZ", "Y vs #theta_{XZ};Y [mm];#theta_{XZ} [#circ]", 500, 0, 256, 150, -10, 10},
                               "Y1", "thetaXZ")};
    auto hZthetaXY {df.Histo2D({"hZthetaXY", "Z vs #theta_{XY};Z [mm];#theta_{XY} [#circ]", 200, 120, 220, 150, -5, 5},
                               "Z1", "thetaXY")};
    auto hZthetaXZ {df.Histo2D({"hZthetaXZ", "Z vs #theta_{XZ};Z [mm];#theta_{XZ} [#circ]", 200, 120, 220, 150, -5, 5},
                               "Z1", "thetaXZ")};
    auto h3d {df.Histo3D({"h3D", "Emittance histogram;Y [mm];#theta_{XY} [#circ];#theta_{XZ} [#circ]", 160, 80, 160,
                          150, -5, 5, 150, -5, 5},
                         "Y1", "thetaXY", "thetaXZ")};

    // Trajectories plot
    double maxxy {257};
    double maxz {320};
    int nbinsxy {200};
    int nbinsz {250};
    ROOT::TThreadedObject<TH2D> hxz {"hxz", "XZ trajectories;X [mm];Z [mm]", nbinsxy, 0, maxxy, nbinsz, 0, maxz};
    ROOT::TThreadedObject<TH2D> hyz {"hyz", "YZ trajectories;Y [mm];Z [mm]", nbinsxy, 0, maxxy, nbinsz, 0, maxz};
    df.Foreach(
        [&](ActRoot::Line& line)
        {
            for(int b = 1; b <= nbinsxy; b++)
            {
                auto x {hxz.Get()->GetXaxis()->GetBinCenter(b)};
                auto pos {line.MoveToX(x)};
                hxz.Get()->Fill(pos.X(), pos.Z());
                hyz.Get()->Fill(pos.Y(), pos.Z());
            }
        },
        {"Line"});

    // Print statistics
    std::cout << "-> Beginning : " << '\n';
    std::cout << "   Mean Y : " << hEntry->GetMean(1) << '\n';
    std::cout << "   FWHM Y : " << hEntry->GetStdDev(1) * 2.35 << '\n';
    std::cout << "   Mean Z : " << hEntry->GetMean(2) << '\n';
    std::cout << "   FWHM Z : " << hEntry->GetStdDev(2) * 2.35 << '\n';
    std::cout << "-> End       : " << '\n';
    std::cout << "   Mean Y : " << hEnd->GetMean(1) << '\n';
    std::cout << "   FWHM Y : " << hEnd->GetStdDev(1) * 2.35 << '\n';
    std::cout << "   Mean Z : " << hEnd->GetMean(2) << '\n';
    std::cout << "   FWHM Z : " << hEnd->GetStdDev(2) * 2.35 << '\n';

    // Fit to get width in Z
    hEntryZ->Fit("gaus", "0QM+");
    auto* fit {hEntryZ->GetFunction("gaus")};
    if(fit)
    {
        auto* text {new TPaveText {0.5, 0.6, 0.7, 0.8, "NDC"}};
        text->SetBorderSize(0);
        text->AddText(TString::Format("#sigma = %.2f mm", fit->GetParameter(2)));
        fit->ResetBit(TF1::kNotDraw);
        hEntryZ->GetListOfFunctions()->Add(text);
    }

    // Plot
    auto c = new TCanvas("c", "Emittance", 1000, 1200);
    c->Divide(2, 3);
    c->cd(1);
    hEntry->DrawClone("colz");
    c->cd(2);
    hEnd->DrawClone("colz");
    c->cd(3);
    hThetaXZ->DrawClone();
    c->cd(4);
    hThetaXY->DrawClone();
    c->cd(5);
    hYThetaXY->DrawClone("colz");
    c->cd(6);
    hYThetaXZ->DrawClone("colz");

    auto* c1 {new TCanvas {"c1", "3D canvas"}};
    h3d->DrawClone();

    // Save objects to file
    TString fOutName = "../Outputs/histos_17F_emittance.root";
    std::cout << "Saving histograms in: " << fOutName << std::endl;
    auto file {std::make_unique<TFile>(fOutName, "recreate")};
    hEntry->Write("hBegin");
    hEnd->Write("hEnd");
    h3d->Write("h3d");
    hxz->Write("hTrajXZ");
    hyz->Write("hTrajYZ");
    hEntryZ->Write("hBeginZ");
}