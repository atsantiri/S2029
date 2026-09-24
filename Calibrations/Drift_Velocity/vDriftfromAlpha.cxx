#include "ActCluster.h"
#include "ActCutsManager.h"
#include "ActDataManager.h"
#include "ActMergerData.h"
#include "ActModularData.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActVoxel.h"

#include <ROOT/RDataFrame.hxx>
#include <random>

#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"

#include <filesystem>

bool LineIntersection(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4,
                      double& ix, double& iy)
{
    // Solving A1x + B1y = C1
    //         A2x + B2y = C2
    //  Using Kramer's method

    double det = (y1 - y2) * (x4 - x3) - (y3 - y4) * (x2 - x1);

    // Parallel (or nearly parallel)
    if(std::abs(det) < 1e-12)
        return false;

    ix = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / det;
    iy = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / det;

    // Check that the intersection lies on both segments
    auto between = [](double a, double b, double c)
    { return c >= std::min(a, b) - 1e-12 && c <= std::max(a, b) + 1e-12; };

    if(between(x1, x2, ix) && between(y1, y2, iy) && between(x3, x4, ix) && between(y3, y4, iy))
        return true;

    return false;
}

struct Source
{
    double x;
    double dx;
    double y;
    double dy;
    TH1D* hIx;
    TH1D* hIy;
};

Source findSource(std::vector<float> fx, std::vector<float> fy, std::vector<float> lx, std::vector<float> ly,
                  int samples = 1e5)
{
    auto hIx = new TH1D {"hIx", "Source Location x [pad]", 150, -50, 0};
    auto hIy = new TH1D {"hIy", "Source Location y [pad]", 150, 20, 60};
    hIx->SetDirectory(nullptr);
    hIy->SetDirectory(nullptr);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, fx.size() - 1);
    for(size_t k = 0; k < samples; ++k)
    {
        size_t i = dist(gen);
        size_t j = dist(gen);
        // avoid comparing a line with itself
        while(j == i)
            j = dist(gen);
        double ix, iy;
        if(LineIntersection(fx[i], fy[i], lx[i], ly[i], fx[j], fy[j], lx[j], ly[j], ix, iy))
        {
            hIx->Fill(ix);
            hIy->Fill(iy);
        }
    }
    TF1* fitIx = new TF1("fitIx", "gaus", -32, -26);
    TF1* fitIy = new TF1("fitIy", "gaus", 35, 45);
    hIx->Fit(fitIx, "Q0R");
    hIy->Fit(fitIy, "Q0R");
    double xSource = fitIx->GetParameter(1);
    double sigmaIx = fitIx->GetParameter(2);
    std::cout << "X intersection mean  = " << xSource << ", sigma = " << sigmaIx << std::endl;
    double ySource = fitIy->GetParameter(1);
    double sigmaIy = fitIy->GetParameter(2);
    std::cout << "Y intersection mean  = " << ySource << ", sigma = " << sigmaIy << std::endl;

    return {xSource, sigmaIx, ySource, sigmaIy, hIx, hIy};
}

std::pair<double, double> vDriftfromAlpha(int run = 62, bool plotting = true)
{
    ROOT::EnableImplicitMT();

    auto d {ROOT::RDataFrame("GETTree", TString::Format("../../RootFiles/Cluster/Clusters_Run_00%d.root", run))};

    // Gate on events with only one cluster
    auto df {d.Filter(
        [](ActRoot::TPCData& data)
        {
            auto size {data.fClusters.size() == 1};
            if(!size)
                return false;
            return data.fClusters.front().GetSizeOfVoxels() >= 30;
        },
        {"TPCData"})};

    // Define last point of cluster in x y z, as the projection of the alpha track
    auto def = df.Define("fLastPoint",
                         [](ActRoot::TPCData& d)
                         {
                             auto cluster {d.fClusters[0]};
                             auto line {cluster.GetRefToLine()};
                             auto dir {line.GetDirection()};
                             cluster.SortAlongDir(dir);
                             auto lastVoxel {cluster.GetRefToVoxels().back()};
                             auto projectionPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
                             return projectionPointLine;
                         },
                         {"TPCData"})
                   .Define("fOtherPoint",
                           [](ActRoot::TPCData& d)
                           {
                               auto cluster {d.fClusters[0]};
                               auto line {cluster.GetRefToLine()};
                               auto otherPoint {line.MoveToX(-50)};
                               return otherPoint;
                           },
                           {"TPCData"})
                   .Define("fLastX", "fLastPoint.X()")
                   .Define("fLastY", "fLastPoint.Y()")
                   .Define("fOtherX", "fOtherPoint.X()")
                   .Define("fOtherY", "fOtherPoint.Y()");

    auto source = findSource(*def.Take<float>("fOtherX"), *def.Take<float>("fOtherY"), *def.Take<float>("fLastX"),
                             *def.Take<float>("fLastY"));

    auto dfDrift =
        def.Define("fDeltaZ",
                   [&](ActRoot::TPCData& d)
                   {
                       if(d.fClusters.size() != 1)
                           return -1000.;
                       else
                       {
                           auto cluster {d.fClusters[0]};
                           auto line {cluster.GetRefToLine()};
                           auto dir {line.GetDirection()};
                           cluster.SortAlongDir(dir);
                           // auto firstVoxel {cluster.GetRefToVoxels().front()};
                           // auto projectionFirstPointLine {line.ProjectionPointOnLine(firstVoxel.GetPosition())};
                           auto lastVoxel {cluster.GetRefToVoxels().back()};
                           auto projectionLastPointLine {line.ProjectionPointOnLine(lastVoxel.GetPosition())};
                           auto zSource {line.MoveToX(source.x).Z()};
                           double deltaZ = projectionLastPointLine.Z() - zSource;
                           return deltaZ * 0.32; // Conversion factor from btb to micro seconds
                       }
                   },
                   {"TPCData"})
            .Define("fLxy",
                    [&](ActRoot::TPCData& d)
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
                            double lxy = TMath::Sqrt(TMath::Power(projectionPointLine.X() - source.x, 2) +
                                                     TMath::Power(projectionPointLine.Y() - source.y, 2));
                            return (lxy * 2) / 10; // Conversion factor from pads to cm
                        }
                    },
                    {"TPCData"})
            .Define("fDeltaZSquare", "fDeltaZ * fDeltaZ")
            .Define("fLxySquare", "fLxy * fLxy");

    auto graphDrift = dfDrift.Graph("fDeltaZ", "fLxy");
    graphDrift->SetTitle("Delta Z vs Lxy;#Deltat [#mus]; #Deltaxy [cm]");

    // Good data cut
    auto cutFile = TString::Format("./Inputs/cut_DriftVelocity_GoodAlphaEvents_%d.root", run);
    if(!std::filesystem::exists(cutFile.Data()))
    {
        TCanvas* c0 = new TCanvas("c0", "Missing cuts", 900, 600);
        graphDrift->GetXaxis()->SetRangeUser(-20, 30);
        graphDrift->GetYaxis()->SetRangeUser(0, 20);
        graphDrift->DrawClone("AP");
        std::cout << "make good events cut and save as " << cutFile << std::endl;
        return {};
    }

    ActRoot::CutsManager<std::string> cuts;
    cuts.ReadCut("goodEvents", cutFile.Data());
    // Filter good events
    auto dff = dfDrift.Filter([&](double lxy, double deltaZ) { return cuts.IsInside("goodEvents", deltaZ, lxy); },
                              {"fLxy", "fDeltaZ"});

    // Peak cuts
    auto peakcutFile = TString::Format("./Inputs/cut_low_%d.root", run);
    if(!std::filesystem::exists(peakcutFile.Data()))
    {
        TCanvas* c0 = new TCanvas("c0", "Missing cuts", 900, 600);
        auto gLinear = dff.Graph("fDeltaZSquare", "fLxySquare");
        gLinear->SetTitle("Delta Z^2 vs Lxy^2;(#Deltat)^{2} [#mus^{2}];(#Deltaxy)^{2} [cm^{2}]");
        gLinear->DrawClone("AP");
        std::cout << "make peak cuts and save as " << peakcutFile << std::endl;
        return {};
    }
    cuts.ReadCut("low", TString::Format("./Inputs/cut_low_%d.root", run).Data());
    cuts.ReadCut("mid", TString::Format("./Inputs/cut_mid_%d.root", run).Data());
    cuts.ReadCut("top", TString::Format("./Inputs/cut_top_%d.root", run).Data());
    auto dfLow = dff.Filter([&](double lxy2, double deltaZ2) { return cuts.IsInside("low", deltaZ2, lxy2); },
                            {"fLxySquare", "fDeltaZSquare"});
    auto dfMid = dff.Filter([&](double lxy2, double deltaZ2) { return cuts.IsInside("mid", deltaZ2, lxy2); },
                            {"fLxySquare", "fDeltaZSquare"});
    auto dfTop = dff.Filter([&](double lxy2, double deltaZ2) { return cuts.IsInside("top", deltaZ2, lxy2); },
                            {"fLxySquare", "fDeltaZSquare"});

    // Make Graphs
    auto fitPeak = [](auto& df, const char* title)
    {
        auto graph = df.Graph("fDeltaZSquare", "fLxySquare");
        graph->SetTitle(title);
        graph->Fit("pol1", "Q");

        return std::make_pair(graph, graph->GetFunction("pol1"));
    };
    auto [glow, f1] = fitPeak(dfLow, "Delta Z^2 vs Lxy^2 (first peak);#Delta Z^2 [#mus^2];Lxy^2 [cm^2]");
    auto [gmid, f2] = fitPeak(dfMid, "Delta Z^2 vs Lxy^2 (second peak);#Delta Z^2 [#mus^2];Lxy^2 [cm^2]");
    auto [gtop, f3] = fitPeak(dfTop, "Delta Z^2 vs Lxy^2 (third peak);#Delta Z^2 [#mus^2];Lxy^2 [cm^2]");

    double v1 = TMath::Sqrt(-f1->GetParameter(1));
    double e1 = TMath::Sqrt(f1->GetParError(1));
    double v2 = TMath::Sqrt(-f2->GetParameter(1));
    double e2 = TMath::Sqrt(f2->GetParError(1));
    double v3 = TMath::Sqrt(-f3->GetParameter(1));
    double e3 = TMath::Sqrt(f3->GetParError(1));

    // Plot is running StandAlone
    if(plotting)
    {
        auto hLast = def.Histo2D({"hLast", "LastPoint XY;X [pads];Y [pads]", 1000, -100, 150, 1000, -100, 150},
                                 "fLastX", "fLastY");
        auto hOther = def.Histo2D({"hOther", "OtherPoint XY;X [pads];Y [pads]", 1000, -100, 150, 1000, -100, 150},
                                  "fOtherX", "fOtherY");
        TCanvas* c = new TCanvas("c", "Points XY / Find Source Location", 1200, 900);
        c->cd();
        TPad* p1 = new TPad("p1", "", 0, 0.5, 1, 1);
        p1->Draw();
        TPad* p2 = new TPad("p2", "", 0, 0, 0.5, 0.5);
        p2->Draw();
        TPad* p3 = new TPad("p3", "", 0.5, 0, 1, 0.5);
        p3->Draw();
        p1->cd();
        hLast->DrawClone("colz");
        hOther->DrawClone("same");
        int counter = 0;
        def.Foreach(
            [&](float otherX, float otherY, float lastX, float lastY)
            {
                counter++;
                auto line = new TLine(otherX, otherY, lastX, lastY);
                line->SetLineColorAlpha(kBlue, 0.3); // transparente para ver cruces
                if(counter % 50 == 0 && lastX > 5 && lastX < 60)
                    line->Draw("same");
            },
            {"fOtherX", "fOtherY", "fLastX", "fLastY"});
        p2->cd();
        source.hIx->Draw();
        p3->cd();
        source.hIy->Draw();

        TCanvas* c1 = new TCanvas("c1", "Delta Z vs Lxy", 900, 600);
        c1->DivideSquare(2);
        c1->cd(1);
        graphDrift->GetXaxis()->SetRangeUser(-20, 30);
        graphDrift->GetYaxis()->SetRangeUser(0, 20);
        graphDrift->DrawClone("AP");
        cuts.DrawCut("goodEvents");

        c1->cd(2);
        auto gLinear = dff.Graph("fDeltaZSquare", "fLxySquare");
        gLinear->SetTitle("Delta Z^2 vs Lxy^2;(#Deltat)^{2} [#mus^{2}];(#Deltaxy)^{2} [cm^{2}]");
        gLinear->DrawClone("AP");

        f1->DrawClone("same");
        cuts.DrawCut("low");
        f2->SetLineColor(kGreen);
        f2->DrawClone("same");
        cuts.DrawCut("mid");
        f3->SetLineColor(kBlue);
        f3->DrawClone("same");
        cuts.DrawCut("top");
        // Text of the fit parameters

        auto t1 = new TLatex(60, 200, TString::Format("First peak: Vdrift = %.2f#pm%.2f ", v1, e1));
        auto t2 = new TLatex(60, 180, TString::Format("Second peak: Vdrift = %.2f#pm%.2f", v2, e2));
        auto t3 = new TLatex(60, 160, TString::Format("Third peak: Vdrift = %.2f#pm%.2f", v3, e3));
        t1->DrawClone();
        t2->DrawClone();
        t3->DrawClone();
    }

    double vdrift = TMath::Mean(3, (double[]) {v1, v2, v3});
    double err = TMath::Sqrt(e1 * e1 + e2 * e2 + e3 * e3) / 3.0;
    return {vdrift, err};
}
