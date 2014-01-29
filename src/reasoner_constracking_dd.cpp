#include <opengm/inference/lpcplex.hxx>

#include "pgmlink/reasoner_constracking_dd.h"

pgmlink::DualDecompositionConservationTracking::DualDecompositionConservationTracking(
        unsigned int max_number_objects,
        boost::function<double (const Traxel&, const size_t)> detection,
        boost::function<double (const Traxel&, const size_t)> division,
        boost::function<double (const double)> transition,
        double forbidden_cost,
        double ep_gap,
        bool with_tracklets,
        bool with_divisions,
        boost::function<double (const Traxel&)> disappearance_cost_fn,
        boost::function<double (const Traxel&)> appearance_cost_fn,
        bool with_misdetections_allowed,
        bool with_appearance,
        bool with_disappearance,
        double transition_parameter,
        bool with_constraints
        )
    : pgmlink::ConservationTracking(
          max_number_objects,
          detection,
          division,
          transition,
          forbidden_cost,
          ep_gap,
          with_tracklets,
          with_divisions,
          disappearance_cost_fn,
          appearance_cost_fn,
          with_misdetections_allowed,
          with_appearance,
          with_disappearance,
          transition_parameter,
          with_constraints),
      hypotheses_graph_(NULL),
      dd_optimizer_(NULL)
{}

pgmlink::DualDecompositionConservationTracking::~DualDecompositionConservationTracking()
{

}

void pgmlink::DualDecompositionConservationTracking::add_constraint(
        pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType* optimizer,
        std::vector<std::size_t>::iterator ids_begin,
        std::vector<std::size_t>::iterator ids_end,
        std::vector<int>::iterator coeffs_begin,
        int lower, int higher, const char* name)
{
    optimizer->addConstraint(ids_begin, ids_end, coeffs_begin, lower, higher, name);
}

void pgmlink::DualDecompositionConservationTracking::configure_hard_constraints(const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType& subGM,
                                                                                pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType& optimizer)
{
    assert(hypotheses_graph_ != NULL);
    //subGM.
    LOG(logDEBUG) << "Found subproblem with " << subGM.numberOfVariables() << " variables";

    pgmlink::ConservationTracking::add_constraints(*hypotheses_graph_, boost::bind(
                                                       &pgmlink::DualDecompositionConservationTracking::add_constraint,
                                                       this, &optimizer, _1, _2, _3, _4, _5, _6));
}

void pgmlink::DualDecompositionConservationTracking::infer()
{
    if (!with_constraints_) {
        opengm::hdf5::save(optimizer_->graphicalModel(), "./conservationTracking.h5", "conservationTracking");
        throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. The conservation tracking factor graph has been saved to file");
    }

    GraphicalModelType* model = pgm_->Model();
    DualDecompositionSubGradient::Parameter dd_parameter;
    dd_parameter.decompositionId_ = DualDecompositionSubGradient::Parameter::TREE;
    dd_parameter.subPara_.verbose_ = true;
    dd_parameter.subPara_.integerConstraint_ = true;
    dd_parameter.subPara_.epGap_ = ep_gap_;

    dd_optimizer_ = new DualDecompositionSubGradient(*model, dd_parameter, boost::bind(&pgmlink::DualDecompositionConservationTracking::configure_hard_constraints, this, _1, _2));

    DualDecompositionSubGradient::VerboseVisitorType visitor;
    opengm::InferenceTermination status = dd_optimizer_->infer(visitor);

    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated abnormally");
    }
}

void pgmlink::DualDecompositionConservationTracking::extractSolution(std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> &solution)
{
    opengm::InferenceTermination status = dd_optimizer_->arg(solution);
    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }
}

void pgmlink::DualDecompositionConservationTracking::add_constraints(const pgmlink::HypothesesGraph &g)
{
    hypotheses_graph_ = &g;
}
