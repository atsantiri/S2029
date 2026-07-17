// Function from https://github.com/loopset/S2008/tree/main/Fits
#ifndef DoubleXS_h
#define DoubleXS_h

#include "ActKinematics.h"
#include "ActSRIM.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TString.h"

#include <utility>
#include <vector>

class DoubleXS
{
public:
    using PairType = std::pair<double, double>;

private:
    std::vector<TH1D*> fProjsThetaCM {};
    std::vector<PairType> fIvsThetaCM {};
    std::vector<TH1D*> fProjsECM {};
    std::vector<PairType> fIvsECM {};
    TH2* fOriginal {};
    TH2* fHist {};
    TH2* fEff {};
    ActPhysics::SRIM* fsrim {};
    TProfile2D* fEffProf {};
    TF1* fFuncOmega {};
    TH2* fThick {};
    ActPhysics::Kinematics* fKin {};
    double fNbeams {};
    double fDensity {};
    TString fIsCM {};

public:
    DoubleXS(TH2* hData, TH2* hEff, ActPhysics::SRIM* srim, double nb, double rho, ActPhysics::Kinematics* kin,
             TString isCM = "Lab");

    void Draw();

    void Project(int cThresch = 50);
    void DrawProjectionsThetaCM(const std::function<void(TH1* h)>& apply = nullptr);
    void DrawProjectionsECM(const std::function<void(TH1* h)>& apply = nullptr);
    TH1D* GetProjectionECM(double thetamin, double thetamax);
    void
    WriteInAzureFormat(int idx, const TString& file, TH1D* pout = nullptr, const std::pair<double, double>& ivs = {-1, -1});
    TH2* GetHist() { return fHist; }

private:
    void ApplyEff();
    void ApplySolidAngle();
    void ApplyThickness();
    void ApplyNormalisation();
    double TransformCMtoLab(double ecm);
};

#endif // !DoubleXS_h
