// I want to test how the Ex would look like if instead of beam energy I was using the recoil nucleus.
// Thiis a function to be called be Pipe3
#include "ActKinematics.h"

#include "TMath.h"

std::pair<double, double> ReconstructExcitationEnergyFromProducts(const ActPhysics::Kinematics& kin, double T3,
                                                                  double theta3LabRads, double T4, double theta4LabRads)
{
    std::pair<double, double> ret;
    const double fm1 {kin.GetParticle(1).GetMass()};
    const double fm2 {kin.GetParticle(2).GetMass()};
    const double fm3 {kin.GetParticle(3).GetMass()};
    const double fm4 {kin.GetParticle(4).GetMass()};

    // check input
    if(!std::isfinite(T3) || !std::isfinite(T4) || !std::isfinite(theta3LabRads) || !std::isfinite(theta4LabRads) ||
       T3 <= 0. || T4 <= 0.)
        return {std::nan(""), std::nan("")};

    auto sin3 {std::sin(theta3LabRads)};
    auto sin4 {std::sin(theta4LabRads)};
    auto cos3 {std::cos(theta3LabRads)};
    auto cos4 {std::cos(theta4LabRads)};

    if(T4 <= 0. || std::abs(sin4) < 1e-9)
        return ret; // recoil along the beam axis or no energy

    auto p3 {TMath::Sqrt(T3 * (T3 + 2.0 * fm3))};

    // // Exact solution. Works only with perfect numbers
    // auto p4 {p3 * sin3 / sin4}; // from momentum conservation along Y
    // auto M4 {(p4 * p4 - T4 * T4) / (2.0 * T4)};
    // auto Ex = M4 - fm4;
    // auto p1 {p3 * cos3 + p4 * cos4};
    // auto E1 {std::sqrt(p1 * p1 + fm1 * fm1)}; // to compare with EBeam

    // exact solution doesn't work, cause Ex = M4 - m4 is a tiny difference of two huge numbers so it’s sensitive to the
    // energy and angle resolution of both particles

    // fixed point iteration instead, https://en.wikipedia.org/wiki/Fixed-point_iteration
    // we calculate the excitation energy as a function of the excitation and iterate until the two converge
    double Ex {0.};
    double E1;
    double E3 {T3 + fm3};
    for(int i = 0; i < 50; i++)
    {
        auto p4 {std::sqrt(T4 * (T4 + 2.0 * (fm4 + Ex)))};
        auto p1 {p3 * cos3 + p4 * cos4};
        E1 = std::sqrt(p1 * p1 + fm1 * fm1);
        double exNew {E1 + fm2 - E3 - T4 - fm4};
        const bool done {std::abs(exNew - Ex) < 1e-9};
        Ex = exNew;
        if(done)
            break;
        if(i > 20)
            std::cout << "Something's off, fixed point iteration can't converge." << std::endl;
    }
    ret = {Ex, E1 - fm1};
    return ret;
}

void testKinematics()
{
    ActPhysics::Kinematics kp {"17F(p,p)17F"};
    auto test {ReconstructExcitationEnergyFromProducts(kp, 1., 20., 0., 15.)};
    std::cout << test.first << " " << test.second << std::endl;
}