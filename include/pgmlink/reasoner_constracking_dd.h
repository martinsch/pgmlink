#ifndef REASONER_CONSTRACKING_DD_H
#define REASONER_CONSTRACKING_DD_H

#include <boost/function.hpp>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include "ext_opengm/dualdecomp_subgradient_hardconstraint.hxx"
#include "reasoner_constracking.h"
#include "pgm.h"

namespace pgmlink {
class Traxel;

class DualDecompositionConservationTracking : public ConservationTracking
{
public:
    typedef pgm::OpengmModelDeprecated::Energy ValueType;
    typedef pgm::OpengmModelDeprecated::ogmGraphicalModel GraphicalModelType;

    typedef opengm::DDDualVariableBlock< marray::Marray<ValueType> > DualBlockType;
    typedef opengm::DualDecompositionBase<pgm::OpengmModelDeprecated::ogmGraphicalModel, DualBlockType>::SubGmType DualDecompositionSubGraphType;
    typedef opengm::LPCplex<DualDecompositionSubGraphType, pgm::OpengmModelDeprecated::ogmAccumulator> InfType;
    typedef opengm::DualDecompositionSubGradientWithHardConstraints<GraphicalModelType,InfType,DualBlockType> DualDecompositionSubGradient;

public:
    DualDecompositionConservationTracking(
            unsigned int max_number_objects,
            boost::function<double (const Traxel&, const size_t)> detection,
            boost::function<double (const Traxel&, const size_t)> division,
            boost::function<double (const double)> transition,
            double forbidden_cost = 0,
            double ep_gap = 0.01,
            bool with_tracklets = false,
            bool with_divisions = true,
            boost::function<double (const Traxel&)> disappearance_cost_fn = ConstantFeature(500.0),
            boost::function<double (const Traxel&)> appearance_cost_fn = ConstantFeature(500.0),
            bool with_misdetections_allowed = true,
            bool with_appearance = true,
            bool with_disappearance = true,
            double transition_parameter = 5,
            bool with_constraints = true
            );

    ~DualDecompositionConservationTracking();

    /**
     * Overwrite the infer method and tell OpenGM to use DualDecomposition for tracking
     */
    virtual void infer();

    void configure_hard_constraints(const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType &subGM,
                                    pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType &optimizer);
protected:
    virtual void extractSolution(std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> &solution);
    DualDecompositionSubGradient* dd_optimizer_;
};

}

#endif // REASONER_CONSTRACKING_DD_H
