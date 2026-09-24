
#include "ActMergerData.h"
#include "ActSilMatrix.h"

namespace S2029
{
auto isFront {[](ActRoot::MergerData& m)
              {
                  if(m.fLight.GetNLayers() == 1)
                      if(m.fLight.GetLayer(0) == "f0")
                          return true;
                  return false;
              }};
auto isLeft {[](ActRoot::MergerData& m)
             {
                 if(m.fLight.GetNLayers() == 1)
                     if(m.fLight.GetLayer(0) == "l0")
                         return true;
                 return false;
             }};
auto isRight {[](ActRoot::MergerData& m)
              {
                  if(m.fLight.GetNLayers() == 1)
                      if(m.fLight.GetLayer(0) == "r0")
                          return true;
                  return false;
              }};
ActPhysics::SilMatrix* GetFrontMatrix()
{
    auto* sm {new ActPhysics::SilMatrix {"f0"}};
    sm->Read("~/ACTAR/S2029/Macros/Outputs/f0_matrix.root");
    return sm;
}

ActPhysics::SilMatrix* GetLeftMatrix()
{
    auto* sm {new ActPhysics::SilMatrix {"l0"}};
    sm->Read("~/ACTAR/S2029/Macros/Outputs/l0_matrix.root");
    return sm;
}

ActPhysics::SilMatrix* GetRightMatrix()
{
    auto* sm {new ActPhysics::SilMatrix {"r0"}};
    sm->Read("~/ACTAR/S2029/Macros/Outputs/r0_matrix.root");
    return sm;
}
} // namespace S2029