#ifndef HistConfig_h
#define HistConfig_h

#include "ROOT/RDF/HistoModels.hxx"

#include "TString.h"

namespace HistConfig
{
using namespace ROOT::RDF;

const TH1DModel Ex {"hEx", TString::Format("Excitation energy;E_{x} [MeV];Counts / %.f keV", (4. - (-2.)) / 200 * 1e3), 200, -2, 4};
const TH1DModel dE {"hdE", "Energy Loss;E_{x} [MeV/u];Counts", 200,0,5.5};
const TH2DModel Ecm_dist {"hEcm_dist","Energy along active area;Distance [mm];E_{CM} [MeV/u]",200,0,256,200,0,4.0};

const TH2DModel Kin {"hKin", "Kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 250, 0, 60, 250, 0, 15};
const TH2DModel KinEl {"hKinEl", "Kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 600, 0, 180, 400, 0, 15};
const TH2DModel KinSimu {"hKin", "Simulation kinematics;#theta_{Lab} [#circ];E_{Vertex} [MeV]", 400, 0, 90, 300, 0, 40};

const TH2DModel KinCM {"hKinCM", "CM kinematics;#theta_{CM} [#circ];E_{Vertex} [MeV]", 400, 0, 180, 400, 0, 15};
const TH1DModel RPx {"hRPx", "RPx;X [mm];Counts", 100, -10, 300};
const TH2DModel RP {"hRP", "RP;X [mm];Y [mm]", 200, -10, 300, 200, -10, 300};

const TH2DModel PID {"hPID", "PID;E_{Sil} [MeV];Q_{ave} [mm^{-1}]", 400, 0, 40, 800, 0, 1000};
const TH2DModel PIDTwo {"hPIDTwo", "PID with two silicons;E_{1} [MeV];E_{0} [MeV]", 800, 0, 40, 800, 0, 40};
const TH2DModel SP {"hSP", "SP;X or Y [mm];Z [mm]", 200, -10, 300, 200, -10, 300};
const TH1DModel TL {"hTL", "Track length; TL [mm]", 300, 0, 600};

const TH1DModel ThetaCM {"hThetaCM", "ThetaCM;#theta_{CM} [#circ]", 600, 0, 180};
const TH2DModel ZThetaZ {"hZThetaZ", "Emittance along Z;Z [mm];#theta_{Z} [#circ]", 600, 0, 270, 600, -10, 10};
const TH2DModel YPhiY {"hYPhiY", "Emittance along Y;Y [mm];#phi_{Y} [#circ]", 600, 0, 270, 600, -10, 10};
const TH2DModel ThetaBeam {"hThetaBeam", "#theta_{Beam} against RP.X;RP.X() [mm];#theta_{Beam} [#circ]", 200, -5, 270, 200, -1, 10};
const TH2DModel ExZ {"hExZ", "E_{x} dependence on SP.Z();SP.Z() [mm];E_{x} [MeV]", 200, -10, 300, 200, -10, 20};
const TH2DModel ExThetaCM {"hExThetaCM", "E_{x} vs #theta_{CM};#theta_{CM} [#circ];E_{x} [MeV]", 400, 0, 180, 200, -5, 10};
const TH2DModel ExThetaLab {"hExThetaLab", "E_{x} vs #theta_{Lab};#theta_{Lab} [#circ];E_{x} [MeV]", 400, 0, 60, 200, -10, 20};
const TH2DModel ExRPx {"hExRPX", "E_{x} vs RP.X;RP.X() [mm];E_{x} [MeV]", 200, -10, 300, 200, -5, 15};
const TH2DModel EBeamRPx {"hEBeamRPX", "E_{Beam} vs RP.X;RP.X() [mm];E_{Beam} [MeV]", 200, -10, 300, 200, 0, 80};

const TH2DModel EcnThetaCM {"hEcnThetaCM", "E_{CN} vs ThetaCM; #theta_CM [#circ]; E_{^{*18}Ne} [MeV]", 720, 0, 180, 200, 4, 9};
const TH2DModel Eff2D {"hEff2D", "2D efficiency; #theta_{CM} [#circ];E_{CM} [MeV]", 720, 0, 180, 200, 0, 5};
const TH2DModel ECMECM {"hECMRes", "E_{CM} resolution;E_{CM}^{nom} [MeV];E_{CM} [MeV]", 500, 0, 5, 500, 0, 5};
const TH2DModel ThetaCMLab {"hThetaCMLab", "CM vs Lab correlations;#theta_{CM} [#circ];#theta_{Lab} [#circ]", 400, 0, 180, 400, 0, 180};
const TH2DModel RPxThetaCM {"hRPxThetaCM", "RP.X vs #theta_{CM} correlations;RP.X [mm];#theta_{CM} [#circ]", 200, 0, 300, 200, 0, 180};
const TH1DModel ECN {"hEcn","Compound Nucleus energy;E_{^{18}Ne} [MeV];Counts ", 150,3.5, 9};
const TH1DModel EBeam {"hEBeam", "Beam energy;E_{Beam} [MeV]", 100, 0, 80};
const TH1DModel ELight {"hELight", "Alpha energy;E_{a} [MeV]", 100, 0, 4};
const TH1DModel EHeavy {"hEHeavy", "Heavy energy;E [MeV]", 100, 0, 70};
const TH2DModel ELightTL {"hELight", "Alpha energy vs Track Length;TL [mm];E_{a} [MeV];",100,0,50, 100, 0, 4 };
const TH1DModel PhiCM {"hPhiCM", "PhiCM;#phi_{CM} [#circ]", 600, 0, 180};
const TH1DModel ECM {"hECM", "E_{CM};E_{CM} [MeV];Counts / 25 keV", 100, 0, 5};
const TH2DModel RPxECM {"hRPxECM", "ECM vs RP.X;RP.X [mm];E_{CM} [MeV]", 200, 0, 300, 150, 0, 15};
const TH2DModel ExESi {"hExESi", "Ex vs E_Si; ESi [MeV]; Ex [MeV]", 150, 0, 15, 100, -5, 5};
const TH2DModel EBeamCompare {"hEBeamCompare", "EBeam vs Reconstructed EBeam; EBeam [MeV]; Reconstructed EBeam [MeV]", 100, 0, 80, 100, 0, 80};


template <typename T>
T ChangeTitle(T model, const TString& title, const TString& label = "");
} // namespace HistConfig

template <typename T>
T HistConfig::ChangeTitle(T model, const TString& title, const TString& label)
{
    auto ret {model};
    if(label.Length() > 0)
        ret.fName = model.fName + label;
    TString old {model.fTitle};
    auto end {old.First(';')};
    TString nt {title + old(end, old.Length() - end)};
    ret.fTitle = nt;
    return ret;
}

#endif 
