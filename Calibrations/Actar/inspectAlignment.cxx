#include "ActDataManager.h"
#include "ActTPCData.h"

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TLine.h"

// I want to get two energy spectra of alpha particles. One extraced from the track lengths, and one from the total
// collected charge along the tracks. I want to compare the resolution of the two methods. Track length resolution
// should in pronciple be better than collected charge (Fig.7. Pancin NIMA 2014)

void inspectAlignment()
{
    ActRoot::DataManager dataman {"../../configs/data.conf", ActRoot::ModeType::EReadTPC};
    dataman.SetRuns(62, 62);
    auto chain {dataman.GetChain()};

    ROOT::RDataFrame d {*chain};

    auto df {d.Filter([](ActRoot::TPCData& tpc) { return (tpc.fClusters.size() == 1); }, {"TPCData"})};

    // found during drift velocity calculation
    float xSource {-29.};
    float ySource {39.};

    auto dfPoints =
        df.Define(
              "lastPoint",
              [](ActRoot::TPCData& tpc)
              {
                  auto cl {tpc.fClusters[0]};              // Get the cluster
                  auto l {cl.GetRefToLine()};              // Get the line that was fitted on the cluster
                  auto dir {l.GetDirection()};             // Get its direction
                  cl.SortAlongDir(dir);                    // sort along the line to find first/last voxel
                  auto voxel {cl.GetRefToVoxels().back()}; // get voxel of interest - here the last
                  auto projection {l.ProjectionPointOnLine(
                      voxel.GetPosition())}; // Get the projection of the position of the first/last voxel on the line
                  return projection;
              },
              {"TPCData"})
            .Define(
                "firstPoint",
                [](ActRoot::TPCData& tpc)
                {
                    auto cl {tpc.fClusters[0]};               // Get the cluster
                    auto l {cl.GetRefToLine()};               // Get the line that was fitted on the cluster
                    auto dir {l.GetDirection()};              // Get its direction
                    cl.SortAlongDir(dir);                     // sort along the line to find first/last voxel
                    auto voxel {cl.GetRefToVoxels().front()}; // get voxel of interest - here the first
                    auto projection {l.ProjectionPointOnLine(
                        voxel.GetPosition())}; // Get the projection of the position of the first/last voxel on the line
                    return projection;
                },
                {"TPCData"})
            .Define("lastPX", "lastPoint.X()")
            .Define("lastPY", "lastPoint.Y()")
            .Define("firstPX", "firstPoint.X()")
            .Define("firstPY", "firstPoint.Y()")
            .Define("Lxy",
                    [&](ROOT::Math::XYZPointF& lastPoint)
                    {
                        double lxy = TMath::Sqrt(TMath::Power(lastPoint.X() - xSource, 2) +
                                                 TMath::Power(lastPoint.Y() - ySource, 2));
                        return lxy * 2; // Conversion factor from pads to mm
                    },
                    {"lastPoint"})
            .Define("dZ",
                    [&](ActRoot::TPCData& tpc, ROOT::Math::XYZPointF& firstPoint, ROOT::Math::XYZPointF& lastPoint)
                    {
                        auto cl {tpc.fClusters[0]};
                        auto l {cl.GetRefToLine()};
                        l.AlignUsingPoint(firstPoint);
                        auto zSource {l.MoveToX(xSource).Z()};
                        double deltaZ = lastPoint.Z() - zSource;
                        return deltaZ * 0.32; // Conversion factor from bin time bucket to micro seconds. 4 time
                                              // buckets in 1 bin time bucket. 1 time bucket is 12.5MHz.
                    },
                    {"TPCData", "firstPoint", "lastPoint"});

    // Plot first and last points to inspect
    auto c = new TCanvas("c", "c", 1000, 1000);
    c->DivideSquare(4);
    c->cd(1);
    auto hLast =
        dfPoints.Histo2D({"hLast", "XY;X [pads];Y [pads]", 1000, -10, 120, 1000, -10, 120}, "lastPX", "lastPY");
    auto hFirst =
        dfPoints.Histo2D({"hFirst", "XY;X [pads];Y [pads]", 1000, -10, 120, 1000, -10, 120}, "firstPX", "firstPY");
    hLast->DrawClone("colz");
    hFirst->DrawClone("same");

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
        {"firstPX", "firstPY", "lastPX", "lastPY"});


    c->cd(2);
    auto gdZLxy = dfPoints.Graph("dZ", "Lxy");
    gdZLxy->SetTitle("Delta Z vs Lxy;#Delta Z [#mus]; Lxy [mm]");
    gdZLxy->GetXaxis()->SetRangeUser(-20, 10);
    gdZLxy->GetYaxis()->SetRangeUser(-0, 200);
    gdZLxy->DrawClone("AP");

    // Get horizontal tracks so I don't need a Z correction and define a total accumulated charge column
    auto dfconstZ = dfPoints.Filter([&](double dZ) { return (dZ >= -2 && dZ < 2); }, {"dZ"})
                            .Define("Qtot",
                            [](ActRoot::TPCData& tpc){
                                auto cl {tpc.fClusters[0]};
                                auto vxs {cl.GetVoxels()};
                                auto qtot {0.};
                                for (const auto &vx: vxs)
                                    qtot += vx.GetCharge();
                                return qtot;
                            },{"TPCData"});


    auto hLxy {dfconstZ.Histo1D({"hLxy", "Alphas for dZ~0 ;Range [mm];Counts", 100, 0, 200}, "Lxy")};
    c->cd(3);
    hLxy->DrawClone();
    auto hQtot {dfconstZ.Histo1D({"hQtot", "Alphas for dZ~0 ;Total Charge [a.u.];Counts", 250, 0, 120000}, "Qtot")};
    c->cd(4);
    hQtot->DrawClone();

    c->SaveAs("alphas_chargeResolution_trackLengthResolution.png");
}