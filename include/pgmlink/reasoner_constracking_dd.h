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

    void configure_hard_constraints(const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType &subGM, size_t sub_gm_index,
                                    pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType &optimizer);
    void add_constraint(pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType *optimizer, size_t sub_gm_index,
                        std::vector<std::size_t>::iterator ids_begin,
                        std::vector<std::size_t>::iterator ids_end,
                        std::vector<int>::iterator coeffs_begin,
                        int lower, int higher, const char *name);
    void debug_graph_output(GraphicalModelType* model);
    void constraint_debug_output(const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType& subGM);
protected:
    virtual void extractSolution(std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> &solution);
    virtual void add_constraints( const HypothesesGraph& );
    virtual size_t cplex_id(size_t opengm_id, size_t state);


    DualDecompositionSubGradient* dd_optimizer_;
    const HypothesesGraph* hypotheses_graph_;
    size_t current_sub_gm_id_;
    DualDecompositionSubGradient::InfType* current_sub_optimizer_;
};

}

#endif // REASONER_CONSTRACKING_DD_H
