#include "ActVAction.h"

namespace ActAlgorithm
{
    class L1Filter : public VAction
    {
    public:
        double fTLmax{}; //!< max TL of noise-like clusters to remove
        double fdZmin{}; //!< minimum dZ distance between noise cluster and RP
        // double fOutOfConeVoxels {}; //!< max number of voxels outside the beam cone to be considered noise of beam spread 
        // double fBeamMaxSpread {}; //!< max spread of beam

    public:
        L1Filter() : VAction("L1Filter") {}

        void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
        void Run() override;
        void Print() const override;

    // private:
    //     struct BeamCone
    //     {
    //         std::vector<std::pair<double, double>> yzWidth;
    //         double x0, xMax, dxStep;
    //         ActRoot::Line line;
    //         bool Contains(ActRoot::Voxel &v) const;
    //         double GetDistfromCone(ActRoot::Voxel &v) const;
    //     };

    //     BeamCone GetBeamCone(const ActRoot::Cluster &c);
    };
} // namespace ActAlgorithm
