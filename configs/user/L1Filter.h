#include "ActVAction.h"

namespace ActAlgorithm
{
    class L1Filter : public VAction
    {
    public:
        double fTLmax{}; //!< max TL of noise-like clusters to remove
        double fdZmin{}; //!< minimum dZ distance between noise cluster and RP

    public:
        L1Filter() : VAction("L1Filter") {}

        void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
        void Run() override;
        void Print() const override;
    };
} // namespace ActAlgorithm
