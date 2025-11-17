#include "ActDataManager.h"
#include "ActLine.h"
#include "ActModularData.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"


void beamEmittance()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EFilter};
    auto chain {dataman.GetChain()};
    auto chain2 {dataman.GetChain(ActRoot::ModeType::EReadSilMod)};
    chain->AddFriend(chain2.get());

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df {*chain};
    // df.Describe().Print();

    auto dff {df
            .Filter([](ActRoot::ModularData& m) { return m.Get("GATCONF") == 8; }, {"ModularData"}) // gate on CFA_DIV
            .Filter("fClusters.fIsBeamLike.size() == 1") // gate on single track events
            .Define("Line",
                    [](ActRoot::TPCData& tpc)
                    {
                        auto cl {tpc.fClusters[0]};
                        auto l {cl.GetRefToLine()};
                        return l;
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

    auto c = new TCanvas("c", "",1000,1200);
    // auto hRange {dff.Histo1D({"hRange","Range in ACTAR;X [mm];Counts",128,0,128},"Range")};
    // hRange->DrawClone();
    c->Divide(2,3);
    c->cd(1);
    auto hEntry {dff.Histo2D({"hEntry", "X = 0 mm;Y [pad];Z [btb]", 128, 0, 128, 128, 0, 128}, "Y1", "Z1")};
    hEntry->DrawClone("colz");
    gPad->SetLogz();
    c->cd(2);
    auto hEnd {dff.Histo2D({"hEnd", "X = 100 mm;Y [pad];Z [btb]", 128, 0, 128, 128, 0, 128}, "Y2", "Z2")};
    hEnd->DrawClone("colz");
    gPad->SetLogz();
    c->cd(3);
    auto hThetaXZ {dff.Histo1D({"hThetaXZ", "XZ Angle;;", 250, -5, 5}, "thetaXZ")};
    hThetaXZ->DrawClone();
    c->cd(4);
    auto hThetaXY {dff.Histo1D({"hThetaXY", "XY Angle;;", 250, -5, 5}, "thetaXY")};
    hThetaXY->DrawClone();
    c->cd(5);
    auto hYThetaXY {dff.Histo2D({"hYThetaXY", "Y vs #theta_{XY};Y [pad];#theta_{XY} [#circ]", 128, 0, 128, 150, -10,
    10},"Y1","thetaXY")}; 
    hYThetaXY->DrawClone("colz"); 
    gPad->SetLogz();
    c->cd(6); 
    auto hYThetaXZ {dff.Histo2D({"hYThetaXZ", "Y vs #theta_{XZ};Y [pad];#theta_{XZ} [#circ]", 128, 0, 128, 150, -10, 10},"Y1","thetaXZ")};
    hYThetaXZ->DrawClone("colz");
    gPad->SetLogz();


    // Print statistics
    std::cout << "-> Beginning : " << '\n';
    std::cout << "   FWHM Y : " << hEntry->GetStdDev(1) * 2.35 << '\n';
    std::cout << "   FWHM Z : " << hEntry->GetStdDev(2) * 2.35 << '\n';
    std::cout << "-> End       : " << '\n';
    std::cout << "   FWHM Y : " << hEnd->GetStdDev(1) * 2.35 << '\n';
    std::cout << "   FWHM Z : " << hEnd->GetStdDev(2) * 2.35 << '\n';
}