// Code from IBC fro s2384 to calculate drift velocity
// expended to look at multiple files with different voltages

#include "ActDataManager.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActModularData.h"
#include "ActModularData.h"
#include "ActMergerData.h"
#include "ActTPCData.h"
#include "ActCluster.h"
#include "ActVoxel.h"
#include "ActCutsManager.h"

#include <ROOT/RDataFrame.hxx>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"

double vDriftfromAlpha(int run)
{
    ROOT::EnableImplicitMT();
    // Get the data run 123
    auto df{ROOT::RDataFrame("GETTree", TString::Format("../../RootFiles/Cluster/Clusters_Run_%04d.root",run))};
    auto vdrift {0.};
    // df.Describe().Print();
    std::cout<<"Processing Run "<<run<<std::endl;

    // Define last point of cluster in x y z, as the projection of the alpha track
    auto dfLastPoint = df.Define("fLastPoint", [](ActRoot::TPCData &d)
                                 {
        if(d.fClusters.size() != 1)
        {
            ROOT::Math::XYZPointF lastPoint {-1000, -1000, -1000};
            return lastPoint;
        }
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            return projectionPointLine;
        } }, {"TPCData"});

    // Define other point
    auto dfOtherPoint = dfLastPoint.Define("fOtherPoint", [](ActRoot::TPCData &d)
    {
        if(d.fClusters.size() != 1)
        {
            ROOT::Math::XYZPointF otherPoint {-1000, -1000, -2000};
            return otherPoint;
        }
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto otherPoint {line.MoveToX(-50)};
            return otherPoint;
        } }, {"TPCData"});

    // Define the coordinates of the points
    auto dfXY = dfOtherPoint
                    .Define("fLastX", "fLastPoint.X()")
                    .Define("fLastY", "fLastPoint.Y()")
                    .Define("fOtherX", "fOtherPoint.X()")
                    .Define("fOtherY", "fOtherPoint.Y()");

    // Plot
    // auto c = new TCanvas("c", "Points XY", 1000, 500);
    // auto hLast = dfXY.Histo2D({"hLast", "LastPoint XY;X [pads];Y [pads]", 1000, -100, 100, 1000, -100, 100}, "fLastX", "fLastY");
    // auto hOther = dfXY.Histo2D({"hOther", "OtherPoint XY;X [pads];Y [pads]", 1000, -100, 100, 1000, -100, 100}, "fOtherX", "fOtherY");
    // hLast->DrawClone("colz");
    // hOther->DrawClone("same");

    // Create the lines and draw them with the foreach
    // int counter = 0;
    // dfXY.Foreach(
    //     [&](float otherX, float otherY, float lastX, float lastY)
    //     {
    //         counter++;
    //         auto line = new TLine(otherX, otherY, lastX, lastY);
    //         line->SetLineColorAlpha(kBlue, 0.3); // transparente para ver cruces
    //         if (counter % 50 == 0 && lastX > 5 && lastX < 60)
    //             line->Draw("same");
    //     },
    //     {"fOtherX", "fOtherY", "fLastX", "fLastY"});

    // With the plot I guess the alpha source is at (-27, 41)
    float xSource{-28.};
    float ySource{39.};

    auto dfDrift = dfXY.Define("fDeltaZ", [&](ActRoot::TPCData &d)
                               {
        if(d.fClusters.size() != 1)
            return -1000.;
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            //auto firstVoxel {cluster.GetRefToVoxels().front()};
            //auto projectionFirstPointLine {line.ProjectionPointOnLine(firstVoxel.GetPosition())};
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionLastPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            auto zSource {line.MoveToX(xSource).Z()};
            double deltaZ = projectionLastPointLine.Z() - zSource;
            return deltaZ *0.32; // Conversion factor from bin time bucket to micro seconds. 4 time buckets in 1 bin time bucket. 1 time bucket is 12.5MHz.
        } }, {"TPCData"})
                       .Define("fLxy", [&](ActRoot::TPCData &d)
                               {
        if(d.fClusters.size() != 1)
            return -1000.;
        else
        {
            auto cluster {d.fClusters[0]};
            auto line {cluster.GetRefToLine()};
            auto dir {line.GetDirection()};
            cluster.SortAlongDir(dir);
            auto lastVoxel {cluster.GetRefToVoxels().back()};
            auto projectionPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
            double lxy = TMath::Sqrt(TMath::Power(projectionPointLine.X() - xSource, 2) + TMath::Power(projectionPointLine.Y() - ySource, 2));
            return lxy * 2; // Conversion factor from pads to mm
        } }, {"TPCData"})
                       .Define("fDeltaZSquare", "fDeltaZ * fDeltaZ")
                       .Define("fLxySquare", "fLxy * fLxy");

    // Plot the DeltaZ and lxy
    auto graphDrift = dfDrift.Graph("fDeltaZ", "fLxy");
    graphDrift->SetTitle("Delta Z vs Lxy;#Delta Z [#mus]; Lxy [mm]");
    graphDrift->GetXaxis()->SetRangeUser(-40,40);
    graphDrift->GetYaxis()->SetRangeUser(-50,300);

    // Linearize the graph
    auto graphDriftLinear = dfDrift.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLinear->SetTitle("Delta Z^2 vs Lxy^2;#Delta Z^2 [#mus^2]; Lxy^2 [mm^2]");


    // auto c1 = new TCanvas("c1", "Delta Z vs Lxy", 1400, 800);
    // c1->DivideSquare(2);
    // c1->cd(1);
    // graphDrift->DrawClone("AP");
    // c1->cd(2);
    // graphDriftLinear->DrawClone("AP");

    // Cuts for good events (no broad region) and for each line
    ActRoot::CutsManager<std::string> cuts;
    // Gas PID
    cuts.ReadCut("goodEvents", "./cut_DriftVelocity_GoodAlphaEvents.root");
    cuts.ReadCut("first", "./cut_firstPeak_run17.root");
    cuts.ReadCut("second", "./cut_secondPeak_run17.root");
    cuts.ReadCut("third", "./cut_thirdPeak_run17.root");

    auto dfFiltered = dfDrift.Filter([&](double lxy, double deltaZ)
                                { return cuts.IsInside("goodEvents", deltaZ, lxy); },
                                {"fLxy", "fDeltaZ"});

    auto graphDriftFiltered = dfFiltered.Graph("fDeltaZ", "fLxy");
    graphDriftFiltered->SetTitle("Delta Z vs Lxy ;#Delta Z [#mus]; Lxy [mm]");
    auto graphDriftFilteredLinear = dfFiltered.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftFilteredLinear->SetTitle("Delta Z^2 vs Lxy^2;#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    // auto c2 = new TCanvas("c2", "Delta Z vs Lxy filtered", 1400, 800);
    // c2->DivideSquare(2);
    // c2->cd(1);
    // graphDriftFiltered->DrawClone("AP");
    // c2->cd(2);
    // graphDriftFilteredLinear->DrawClone("AP");

    // Do graphs for each peak
    auto dfFirst = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("first", deltaZ2, lxy2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineFirst = dfFirst.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineFirst->SetTitle("Delta Z^2 vs Lxy^2 (first peak);#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    graphDriftLineFirst->Fit("pol1");
    auto f1 {graphDriftLineFirst->GetFunction("pol1")};
    f1->SetLineColor(kRed);
    auto dfSecond = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("second", deltaZ2,lxy2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineSecond = dfSecond.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineSecond->SetTitle("Delta Z^2 vs Lxy^2 (second peak);#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    graphDriftLineSecond->Fit("pol1");
    auto f2 {graphDriftLineSecond->GetFunction("pol1")};
    auto dfThird = dfFiltered.Filter([&](double lxy2, double deltaZ2)
                                { return cuts.IsInside("third", deltaZ2,lxy2); },
                                {"fLxySquare", "fDeltaZSquare"});
    auto graphDriftLineThird = dfThird.Graph("fDeltaZSquare", "fLxySquare");
    graphDriftLineThird->SetTitle("Delta Z^2 vs Lxy^2 (third peak);#Delta Z^2 [#mus^2];Lxy^2 [mm^2]");
    graphDriftLineThird->Fit("pol1");
    auto f3 {graphDriftLineThird->GetFunction("pol1")};
    
    // auto c3 = new TCanvas("c3", "Delta Z vs Lxy lines", 2100, 700);
    // c3->DivideSquare(3);
    // c3->cd(1);
    // graphDriftLineFirst->DrawClone("AP");
    // f1->Draw("same");
    // c3->cd(2);
    // graphDriftLineSecond->DrawClone("AP");
    // f2->Draw("same");
    // c3->cd(3);
    // graphDriftLineThird->DrawClone("AP");
    // f3->Draw("same");

    // Draw them also in the filtered plot
    // c2->cd(2);
    // f1->DrawClone("same");
    // f2->SetLineColor(kGreen);
    // f2->DrawClone("same");
    // f3->SetLineColor(kBlue);
    // f3->DrawClone("same");
    // // Text of the fit parameters
    // auto t1 = new TLatex(50, 30000, TString::Format("First peak: Vdrift = %.2f#pm%.2f ", TMath::Sqrt(-f1->GetParameter(1)),TMath::Sqrt(f1->GetParError(1))));
    // auto t2 = new TLatex(50, 28000, TString::Format("Second peak: Vdrift = %.2f#pm%.2f", TMath::Sqrt(-f2->GetParameter(1)),TMath::Sqrt(f2->GetParError(1))));
    // auto t3 = new TLatex(50, 25000, TString::Format("Third peak: Vdrift = %.2f#pm%.2f", TMath::Sqrt(-f3->GetParameter(1)),TMath::Sqrt(f3->GetParError(1))));
    // t1->DrawClone();
    // t2->DrawClone();
    // t3->DrawClone();

    vdrift = TMath::Mean(3,(double[]){TMath::Sqrt(-f1->GetParameter(1)), TMath::Sqrt(-f2->GetParameter(1)),TMath::Sqrt(-f3->GetParameter(1))});

    return vdrift;

}

void multiRuns_vDriftfromAlpha(){
/* List of runs for drift measurements
run 	Vdt (V) 	Vdb (V) 	Vm (V)  deltaV    comment
2 	    5500 	    430 	    430 	5070        GET gain = 120  fC	 
3 	    5500 	    450 	    450 	5050        GET gain = 120  fC	 
4 	    5500 	    450 	    450 	5050        GET gain = 240  fC	 
5 	    5500 	    470 	    470 	5030        GET gain = 240  fC	 
6 	    5500 	    490 	    490 	5010        GET gain = 240  fC	junk - lost connection to one of the CoBos
7 	    5500 	    490 	    490 	5010        GET gain = 240  fC	  	 
8 	    5500 	    490 	    490 	5010        GET gain = 1000 fC 	 
9 	    5500 	    510 	    510 	4990        GET gain = 1000 fC 	 
10 	    5500 	    530 	    530 	4970        GET gain = 1000 fC 	 
11 	    5500 	    550 	    550 	4950        GET gain = 1000 fC 	 
17 	    5500 	    430 	    430 	5070         
18 	    5800 	    430 	    430 	5370         
19 	    6100 	    430 	    430 	5670        1 CoBo failed during the run
20 	    6100 	    430 	    430 	5670         
21 	    5950 	    430 	    430 	5520         
*/
    
    std::map<int, double> data ={
        // { 2, 5070},
        // { 3, 5050},
        // { 4, 5050},
        // { 5, 5030},
        // // { 6, 5010},
        // { 7, 5010},
        // { 8, 5010},
        // { 9, 4990},
        // {10, 4970},
        // {11, 4950},
        {17, 5070},
        {18, 5370},
        // {19, 5670},
        {20, 5670},
        {21, 5520}
    };

    std::vector<double> drifts;
    auto* hDrifts = new TGraph();
    auto nPoints {0};

    for (const auto& m : data){
        auto drift = vDriftfromAlpha(m.first);
        drifts.push_back(drift);
        std::cout<<m.second<<" "<<drift<<std::endl;
        hDrifts->SetPoint(nPoints,m.second,drift);    
        nPoints++;
    }

    auto c3 = new TCanvas("c3", "Drifts vs deltaV", 1400, 800);
    c3->cd();

    hDrifts->SetMarkerStyle(20);
    hDrifts->SetMarkerSize(1.2);
    hDrifts->SetMarkerColor(kBlue);
    hDrifts->SetTitle("Drift vs deltaV;deltaV [V];Drift [mm/us]");

    hDrifts->Draw("AP");
    c3->Update();
    c3->SaveAs("vdrift_runs17-21.png");


}