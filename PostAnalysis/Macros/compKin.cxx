#include "ActKinematics.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
void compKin()
{
    ActPhysics::Kinematics kf {"17F(p,p)@50"};
    ActPhysics::Kinematics ko {"17O(p,p)@65"};

    auto* mg {new TMultiGraph};
    mg->SetTitle("Kin comparison;#theta_{Lab} [#circ];E_{Lab} [MeV]");
    auto* gf {kf.GetKinematicLine3()};
    gf->SetTitle("17F");
    auto* go {ko.GetKinematicLine3()};
    go->SetTitle("17O");
    mg->Add(gf);
    mg->Add(go);

    auto* c0 {new TCanvas {"c0", "Comp kin"}};
    mg->Draw("apl plc pmc");
}
