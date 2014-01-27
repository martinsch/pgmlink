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
          with_constraints)
{}

pgmlink::DualDecompositionConservationTracking::~DualDecompositionConservationTracking()
{

}

void pgmlink::DualDecompositionConservationTracking::infer()
{
    if (!with_constraints_) {
        opengm::hdf5::save(optimizer_->graphicalModel(), "./conservationTracking.h5", "conservationTracking");
        throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. The conservation tracking factor graph has been saved to file");
    }

    GraphicalModelType* model = pgm_->Model();
    dd_optimizer_ = new DualDecompositionSubGradient(*model);

    opengm::InferenceTermination status = dd_optimizer_->infer();

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
