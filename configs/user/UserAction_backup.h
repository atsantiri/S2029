#include "ActVAction.h"

namespace ActAlgorithm
{
class UserAction : public VAction
{
public:
    double fMinLength {};      //!< Min length in X of the cluster set as reference
    double fMaxAngle {};       //!< Max angle for reference cluster in DEGREES
    double fCylinderR {};      //!< Radius of cylinder to compare if clusters overlap

public:
    UserAction() : VAction("UserAction") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm
