#include "ActDataManager.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"

#include "ROOT/RDataFrame.hxx"
#include <random>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"

// I want to fine tune my SRIM tables so the energy of the alphas gets reproduced. Find source location and plot TL vs
// expected E for drift run and from SRIM for different pressures. Aurora's Fig 5.21

ROOT::Math::XYZPointF
findStartVoxel(ActRoot::Cluster& c, bool returnStart) // compute the distance to the center mass of the charge. The one
                                                      // closer to the center mass is the end point
{
    float qtot = 0;
    float xcm = 0., ycm = 0., zcm = 0.;

    for(const auto& v : c.GetVoxels())
    {
        ROOT::Math::XYZPointF pos = v.GetPosition();
        auto q {v.GetCharge()};

        qtot += q;
        xcm += q * pos.X();
        ycm += q * pos.Y();
        zcm += q * pos.Z();
    }

    if(qtot == 0)
        return {};
    xcm /= qtot;
    ycm /= qtot;
    zcm /= qtot;
    ROOT::Math::XYZPointF rcm(xcm, ycm, zcm);
    ROOT::Math::XYZPointF A = c.GetRefToVoxels().front().GetPosition();
    ROOT::Math::XYZPointF B = c.GetRefToVoxels().back().GetPosition();
    float dA = (A - rcm).R();
    float dB = (B - rcm).R();

    ROOT::Math::XYZPointF start = (dA > dB) ? A : B;
    ROOT::Math::XYZPointF end = (dA > dB) ? B : A;

    return returnStart ? start : end;
}

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
// Scale, ScalePoint and calcTl functions from ActMergerDetector and ActLine classes
void ActRoot::Line::Scale(float xy, float z)
{
    // Point
    fPoint.SetX(fPoint.X() * xy);
    fPoint.SetY(fPoint.Y() * xy);
    fPoint.SetZ(fPoint.Z() * z);
    // Direction
    fDirection.SetX(fDirection.X() * xy);
    fDirection.SetY(fDirection.Y() * xy);
    fDirection.SetZ(fDirection.Z() * z);
}

void ScalePoint(ROOT::Math::XYZPointF& point, float xy, float z)
{
    point += ROOT::Math::XYZVector {0.5, 0.5, 0.5};
    point.SetX(point.X() * xy);
    point.SetY(point.Y() * xy);
    point.SetZ(point.Z() * z);
}

double calcTLfromVoxel(ROOT::Math::XYZPointF A, ROOT::Math::XYZPointF B, ActRoot::Line line, double drift)
{
    ActRoot::TPCParameters params;
    double xy = params.GetPadSide();

    ScalePoint(A, xy, drift);
    ScalePoint(B, xy, drift);
    line.Scale(xy, drift);

    auto projBegin {line.ProjectionPointOnLine(A)};
    auto projEnd {line.ProjectionPointOnLine(B)};
    return (projBegin - projEnd).R();
}

void fineTuneSRIMfromAlphas()
{
    ActRoot::InputParser parserDet {"../../configs/detector.conf"};
    auto bl1 {parserDet.GetBlock("Merger")};
    auto drift {bl1->GetDouble("DriftFactor")};

    auto* srim {new ActPhysics::SRIM()};

    auto d {ROOT::RDataFrame("GETTree", "../../RootFiles/Cluster/Clusters_Run_0062.root")};
    auto df {d.Filter([](ActRoot::TPCData& tpc) { return (tpc.fClusters.size() == 1); }, {"TPCData"})};

    auto dfPoints = df.Define("lastPoint",
                              [](ActRoot::TPCData& tpc)
                              {
                                  auto cl {tpc.fClusters[0]};  // Get the cluster
                                  auto l {cl.GetRefToLine()};  // Get the line that was fitted on the cluster
                                  auto dir {l.GetDirection()}; // Get its direction
                                  cl.SortAlongDir(dir);        // sort along the line to find first/last voxel
                                  auto vPos = findStartVoxel(cl, false);
                                  // Get the projection of the position of the first/last voxel on the line
                                  auto projection {l.ProjectionPointOnLine(vPos)};
                                  return projection;
                              },
                              {"TPCData"})
                        .Define("firstPoint",
                                [](ActRoot::TPCData& tpc)
                                {
                                    auto cl {tpc.fClusters[0]};  // Get the cluster
                                    auto l {cl.GetRefToLine()};  // Get the line that was fitted on the cluster
                                    auto dir {l.GetDirection()}; // Get its direction
                                    cl.SortAlongDir(dir);        // sort along the line to find first/last voxel
                                    auto vPos = findStartVoxel(cl, true);
                                    // Get the projection of the position of the first/last voxel on the line
                                    auto projection {l.ProjectionPointOnLine(vPos)};
                                    return projection;
                                },
                                {"TPCData"})
                        .Define("otherPoint",
                                [](ActRoot::TPCData& tpc)
                                {
                                    auto cl {tpc.fClusters[0]};       // Get the cluster
                                    auto l {cl.GetRefToLine()};       // Get the line that was fitted on the cluster
                                    auto otherPoint {l.MoveToX(-50)}; //
                                    return otherPoint;
                                },
                                {"TPCData"})
                        .Define("lastPX", "lastPoint.X()")
                        .Define("lastPY", "lastPoint.Y()")
                        .Define("firstPX", "firstPoint.X()")
                        .Define("firstPY", "firstPoint.Y()")
                        .Define("fOtherX", "otherPoint.X()")
                        .Define("fOtherY", "otherPoint.Y()");

    // Expand lines and find intersection point.
    TCanvas* c = new TCanvas("c", "Find Source Location", 1200, 900);
    TPad* p1 = new TPad("p1", "", 0, 0.5, 1, 1);
    p1->Draw();
    TPad* p2 = new TPad("p2", "", 0, 0, 0.5, 0.5);
    p2->Draw();
    TPad* p3 = new TPad("p3", "", 0.5, 0, 1, 0.5);
    p3->Draw();
    p1->cd();
    auto hLast =
        dfPoints.Histo2D({"hLast", "XY;X [pads];Y [pads]", 1000, -100, 120, 1000, -10, 120}, "lastPX", "lastPY");
    hLast->DrawClone("colz");
    auto hOther =
        dfPoints.Histo2D({"hOther", "XY;X [pads];Y [pads]", 1000, -100, 120, 1000, -10, 120}, "fOtherX", "fOtherY");
    hOther->DrawClone("same");
    // Create the lines and draw them with the foreach
    int counter = 0;
    dfPoints.Foreach(
        [&](float fx, float fy, float lx, float ly)
        {
            counter++;
            auto line = new TLine(fx, fy, lx, ly);
            line->SetLineColorAlpha(kBlue, 0.3);
            if(counter % 50 == 0 && lx > 5 &&
               lx < 60) // counter to downsample how many are drawn and x limits to remove additional noise
                line->Draw("same");
        },
        {"fOtherX", "fOtherY", "lastPX", "lastPY"});
    // Find intersection on a smaple of random lines (checking all of them is too many)
    auto fx = *dfPoints.Take<float>("fOtherX");
    auto fy = *dfPoints.Take<float>("fOtherY");
    auto lx = *dfPoints.Take<float>("lastPX");
    auto ly = *dfPoints.Take<float>("lastPY");
    auto hIx = new TH1D {"hIx", "Source Location x [pad]", 150, -50, 0};
    auto hIy = new TH1D {"hIy", "Source Location y [pad]", 150, 20, 60};

    const int samples = 100000;
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
            auto marker = new TMarker(ix, iy, 1);
            marker->SetMarkerColor(kRed);
            marker->Draw("same");
            hIx->Fill(ix);
            hIy->Fill(iy);
        }
    }
    p2->cd();
    hIx->Draw();
    TF1* fitIx = new TF1("fitIx", "gaus", -35, -25);
    hIx->Fit(fitIx, "QR");
    double sourceX = fitIx->GetParameter(1);
    double sigmaIx = fitIx->GetParameter(2);
    std::cout << "X intersection mean  = " << sourceX << ", sigma = " << sigmaIx << std::endl;

    p3->cd();
    hIy->Draw();
    TF1* fitIy = new TF1("fitIy", "gaus", 35, 45);
    hIy->Fit(fitIy, "QR");
    double sourceY = fitIy->GetParameter(1);
    double sigmaIy = fitIy->GetParameter(2);
    std::cout << "Y intersection mean  = " << sourceY << ", sigma = " << sigmaIy << std::endl;


    // with a better constrained source location, calculate TL from last voxel
    auto dfFinal {dfPoints.Define("TL",
                    [&](ActRoot::TPCData& tpc, ROOT::Math::XYZPointF& lastPoint)
                    {
                        auto cl {tpc.fClusters[0]};
                        auto l {cl.GetRefToLine()};
                        auto sourceZ {l.MoveToX(sourceX).Z()};
                        ROOT::Math::XYZPointF source(sourceX, sourceY, sourceZ);
                        double TL = calcTLfromVoxel(source, lastPoint, l, drift);
                        return TL;
                    },
                                  {"TPCData", "lastPoint"})};

    std::vector<double> energies {{5156.59, 5485.56, 5804.77}};

    TCanvas* c2 = new TCanvas("c2", "TL Dist", 800, 600);
    auto hTL {dfFinal.Histo1D({"hTL", "TL;TL [mm];Counts", 150, 50, 200}, "TL")};
    hTL->DrawClone();


    TF1* fTL = new TF1("fTL", "gaus(0)+gaus(3)+gaus(6)", 120, 180);
    fTL->SetParameters(1300, 140, 2, 1000, 155, 2, 720, 170, 2);
    hTL->Fit(fTL, "RQ");
    hTL->DrawClone();
    std::vector<double> TLvalues;
    for(double TL = 50; TL <= 200; TL += 1)
        TLvalues.push_back(TL);

    std::vector<double> TLs = {fTL->GetParameter(1), fTL->GetParameter(4), fTL->GetParameter(7)};
    std::vector<double> TLerr = {fTL->GetParameter(2), fTL->GetParameter(5), fTL->GetParameter(8)};
    auto c3 = new TCanvas("c3", "Energy vs TL", 1000, 700);
    c3->cd();
    auto* grTL = new TGraphErrors();
    for(size_t i = 0; i < TLs.size(); ++i)
    {
        grTL->SetPoint(i, TLs[i], energies[i]);
        grTL->SetPointError(i, TLerr[i], 0);
    }
    grTL->SetMarkerStyle(20);
    grTL->SetTitle(";TL [mm];Energy [keV]");
    grTL->Draw("AP");
    int colors[] = {kBlue,      kOrange + 7, kGreen + 2,  kRed + 1,  kViolet + 1, kCyan + 1, kMagenta + 1,
                    kAzure + 2, kOrange + 1, kSpring + 5, kPink + 5, kTeal + 3,   kGray + 2};
    int k = 0;
    auto* legend = new TLegend(0.7, 0.15, 0.9, 0.5);
    for(int p = 740; p <= 780; p += 5)
    {
        std::string file = Form("../../Simulation/SRIM/4He_H2-iC4H10_95-5_%dmbar.txt", p);
        auto* srim = new ActPhysics::SRIM();
        srim->ReadTable(Form("HeInGas%d", p), file);

        auto* gr = new TGraph();

        for(size_t i = 0; i < TLvalues.size(); ++i)
        {
            double E = srim->EvalInitialEnergy(Form("HeInGas%d", p), 0, TLvalues[i]) * 1000.;
            gr->SetPoint(i, TLvalues[i], E);
        }
        gr->SetTitle(Form("%s;TL [mm];Energy [keV]", file.c_str()));
        gr->SetLineWidth(2);
        gr->SetLineColor(colors[k++]);
        gr->Draw("L same");
        legend->AddEntry(gr, Form("%d mbar", p), "l");
    }
    legend->Draw();
}