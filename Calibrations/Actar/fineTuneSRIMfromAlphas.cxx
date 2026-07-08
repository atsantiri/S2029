#include "ActDataManager.h"
#include "ActSRIM.h"
#include "ActTPCData.h"
#include "ActTPCParameters.h"

#include "ROOT/RDataFrame.hxx"
#include <random>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMarker.h"

// I want to fine tune my SRIM tables so the energy of the alphas gets reproduced.

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
    srim->ReadTable("HeInGas", "../../Simulation/SRIM/4He_H2-iC4H10_95-5_755mbar.txt");

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
    TF1* fitIx = new TF1("fitIx", "gaus", -32, -26);
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
    auto dfFinal {
        dfPoints
            .Define("TL",
                    [&](ActRoot::TPCData& tpc, ROOT::Math::XYZPointF& lastPoint)
                    {
                        auto cl {tpc.fClusters[0]};
                        auto l {cl.GetRefToLine()};
                        auto sourceZ {l.MoveToX(sourceX).Z()};
                        ROOT::Math::XYZPointF source(sourceX, sourceY, sourceZ);
                        double TL = calcTLfromVoxel(source, lastPoint, l, drift);
                        return TL;
                    },
                    {"TPCData", "lastPoint"})
            .Define("Ene", [&](double TL) { return srim->EvalInitialEnergy("HeInGas", 0, TL) * 1000.; }, {"TL"})
            .Define("TLmin",
                    [&](ActRoot::TPCData& tpc, ROOT::Math::XYZPointF& lastPoint)
                    {
                        auto cl {tpc.fClusters[0]};
                        auto l {cl.GetRefToLine()};
                        auto sourceZ {l.MoveToX(sourceX).Z()};
                        ROOT::Math::XYZPointF source(sourceX+sigmaIx, sourceY, sourceZ);
                        double TL = calcTLfromVoxel(source, lastPoint, l, drift);
                        return TL;
                    },
                    {"TPCData", "lastPoint"})
            .Define("Enemin", [&](double TL) { return srim->EvalInitialEnergy("HeInGas", 0, TL) * 1000.; }, {"TLmin"})
            .Define("TLmax",
                    [&](ActRoot::TPCData& tpc, ROOT::Math::XYZPointF& lastPoint)
                    {
                        auto cl {tpc.fClusters[0]};
                        auto l {cl.GetRefToLine()};
                        auto sourceZ {l.MoveToX(sourceX).Z()};
                        ROOT::Math::XYZPointF source(sourceX-sigmaIx, sourceY, sourceZ);
                        double TL = calcTLfromVoxel(source, lastPoint, l, drift);
                        return TL;
                    },
                    {"TPCData", "lastPoint"})
            .Define("Enemax", [&](double TL) { return srim->EvalInitialEnergy("HeInGas", 0, TL) * 1000.; }, {"TLmax"})};
    TCanvas* c1 = new TCanvas("c1", "Energy Reconstruction", 1200, 900);
    auto hEne {dfFinal.Histo1D({"hEne", "Reconstructed alpha energy; Energy [keV];Counts", 150, 3000, 7000}, "Ene")};
    auto hEneMin {dfFinal.Histo1D({"hEneMin", "Reconstructed alpha energy; Energy [keV];Counts", 150, 3000, 7000}, "Enemin")};
    auto hEneMax {dfFinal.Histo1D({"hEneMax", "Reconstructed alpha energy; Energy [keV];Counts", 150, 3000, 7000}, "Enemax")};
    hEneMax->DrawClone();
    // hEneMin->SetLineColor(kMagenta);
    // hEneMin->DrawClone("same");
    // hEneMax->SetLineColor(kOrange);
    // hEneMax->DrawClone("same");

    std::vector<double> energies {{5156.59, 5485.56, 5804.77}};
    int i = 0;
    for(auto& s : energies)
    {
        TLine* st = new TLine(s, 0, s, 1600 - i * 200);
        st->SetLineColor(2);
        st->SetLineStyle(2);
        st->Draw("same");
        auto t = new TLatex(s + 50, 1450 - i * 200, Form("%.2f", s));
        t->SetTextColor(2);
        t->SetTextSize(0.03);
        t->Draw("same");
        i++;
    }


    // c->cd(2);
    // auto gdZLxy = dfPoints.Graph("dZ", "Lxy");
    // gdZLxy->SetTitle("Delta Z vs Lxy;#Delta Z [#mus]; Lxy [mm]");
    // gdZLxy->GetXaxis()->SetRangeUser(-20, 10);
    // gdZLxy->GetYaxis()->SetRangeUser(-0, 200);
    // gdZLxy->DrawClone("AP");


    // // Get horizontal tracks so I don't need a Z correction and define a total accumulated charge column
    // auto dfconstZ = dfPoints
    //                     .Filter([&](double dZ, float dY)
    //                             { return (dZ >= dzmin && dZ <= dzmax && dY >= dymin && dY <= dymax); }, {"dZ",
    //                             "dY"})
    //                     .Define("Qtot",
    //                             [](ActRoot::TPCData& tpc)
    //                             {
    //                                 auto cl {tpc.fClusters[0]};
    //                                 auto vxs {cl.GetVoxels()};
    //                                 auto qtot {0.};
    //                                 for(const auto& vx : vxs)
    //                                     qtot += vx.GetCharge();
    //                                 return qtot;
    //                             },
    //                             {"TPCData"});


    // auto hLxy {dfconstZ.Histo1D(
    //     {"hLxy", TString::Format("dZ in (%.1f,%.1f) & dY in (%.1f,%.1f);Range [mm];Counts", dzmin, dzmax, dymin,
    //     dymax),
    //      100, 0, 200},
    //     "Lxy")};
    // c->cd(3);

    // TF1* fLxy = new TF1("fLxy", "gaus(0)+gaus(3)+gaus(6)", 120, 180);
    // fLxy->SetParameters(160, 140, 2, 120, 150, 2, 50, 170, 2);
    // hLxy->Fit(fLxy, "Q");

    // hLxy->DrawClone();
    // // if(fLxy)
    // //     fLxy->Draw("same");

    // auto hQtot {dfconstZ.Histo1D({"hQtot",
    //                               TString::Format("dZ in (%.1f,%.1f) & dY in (%.1f,%.1f);Total Charge
    //                               [a.u.];Counts",
    //                                               dzmin, dzmax, dymin, dymax),
    //                               200, 0, 120000},
    //                              "Qtot")};
    // c->cd(4);
    // TF1* fQ = new TF1("fQ", "gaus(0)+gaus(3)+gaus(6)", 65000, 95000);
    // fQ->SetParameters(70, 72000, 1200, 40, 82000, 1400, 25, 90000, 1400);
    // hQtot->Fit(fQ, "Q");

    // hQtot->DrawClone();
    // // if(fQ)
    // //     fQ->Draw("same");

    // c->SaveAs("./Outputs/alphas_chargeResolution_trackLengthResolution.png");

    // std::cout << "Range spectrum resolution: " << std::setprecision(2)
    //           << 2.355 * fLxy->GetParameter(2) / fLxy->GetParameter(1) * 100 << " %, "
    //           << 2.355 * fLxy->GetParameter(5) / fLxy->GetParameter(4) * 100 << " %, and "
    //           << 2.355 * fLxy->GetParameter(8) / fLxy->GetParameter(7) * 100 << " %" << std::endl;


    // std::cout << "Charge spectrum resolution: " << std::setprecision(2)
    //           << 2.355 * fQ->GetParameter(2) / fQ->GetParameter(1) * 100 << " %, "
    //           << 2.355 * fQ->GetParameter(5) / fQ->GetParameter(4) * 100 << " %, and "
    //           << 2.355 * fQ->GetParameter(8) / fQ->GetParameter(7) * 100 << " %" << std::endl;

    // // Express track length as energy with SRIM
    // auto* srim {new ActPhysics::SRIM()};
    // srim->ReadTable("alphasInGas", "../../Simulation/SRIM/4He_H2-iC4H10_95-5_755mbar.txt");

    // auto dfEne = dfconstZ.Define("Ene",
    //                              [&](double Lxy)
    //                              {
    //                                  double eneMeV {srim->EvalInitialEnergy("alphasInGas", 0, Lxy)};
    //                                  return eneMeV * 1000.;
    //                              },
    //                              {"Lxy"});
    // auto hEne {dfEne.Histo1D({"hEne", "Reconstructed alpha energy; Energy [keV];Counts", 150, 3000, 7000},
    // "Ene")}; auto c1 = new TCanvas("c1", "c1"); c1->cd(); TF1* fEne = new TF1("fEne", "gaus(0)+gaus(3)+gaus(6)",
    // 4800, 6000); fEne->SetParameters(100, 5100, 10, 60, 5400, 10, 50, 5700, 10); hEne->Fit(fEne, "Q");
    // hEne->DrawClone();
    // std::cout << "Alpha Energy: " << std::fixed << std::setprecision(0) << fEne->GetParameter(1) << " ("
    //           << std::setprecision(2) << 2.355 * fEne->GetParameter(2) / fEne->GetParameter(1) * 100 << " %), "
    //           << std::setprecision(0) << fEne->GetParameter(4) << " (" << std::setprecision(2)
    //           << 2.355 * fEne->GetParameter(5) / fEne->GetParameter(4) * 100 << " %), and " <<
    //           std::setprecision(0)
    //           << fEne->GetParameter(7) << " (" << std::setprecision(2)
    //           << 2.355 * fEne->GetParameter(8) / fEne->GetParameter(7) * 100 << " %)" << std::endl;
}