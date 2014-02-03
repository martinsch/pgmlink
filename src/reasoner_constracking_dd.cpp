#include <opengm/inference/lpcplex.hxx>
#include <opengm/graphicalmodel/graphicalmodel_factor.hxx>
#include <lemon/graph_to_eps.h>
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
        size_t sub_gm_index,
        std::vector<std::size_t>::iterator ids_begin,
        std::vector<std::size_t>::iterator ids_end,
        std::vector<int>::iterator coeffs_begin,
        int lower, int higher, const char* name)
{
    for(std::vector<std::size_t>::iterator it = ids_begin; it != ids_end; ++it)
    {
        if(*it == (size_t)-1)
        {
            LOG(pgmlink::logDEBUG4) << "Discarding constraint: " << name;
            return;
        }
    }

    optimizer->addConstraint(ids_begin, ids_end, coeffs_begin, lower, higher, name);
}

void pgmlink::DualDecompositionConservationTracking::constraint_debug_output(
        const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType& subGM)
{
    lemon::InDegMap<HypothesesGraph> degree(*hypotheses_graph_);

    // compare subGM to hypotheses_graph
    for(size_t i = 0; i < std::max(subGM.numberOfVariables(), pgm_->Model()->numberOfVariables()); i++)
    {
        if(i < subGM.numberOfVariables())
        {
            size_t num_factors_opengm = subGM.numberOfFactors(i);
            std::cout << "\tVariable " << i << " has " << num_factors_opengm << " factors in OpenGM->Decomposition" << std::endl;
        }

        if(i < pgm_->Model()->numberOfVariables())
        {
            size_t num_factors_pgm = pgm_->Model()->numberOfFactors(i);
            std::cout << "\tVariable " << i << " has " << num_factors_pgm << " factors in PGM" << std::endl;
        }
    }

    for(size_t i = 0; i < std::max(subGM.numberOfFactors(), pgm_->Model()->numberOfFactors()); i++)
    {
        if(i < subGM.numberOfFactors())
        {
            std::cout << "\tFactor " << i << " has " << subGM[i].numberOfVariables() << " variables in OpenGM->Decomposition\n\t\t";

            for(size_t v = 0; v < subGM[i].numberOfVariables(); v++)
            {
                std::cout << subGM[i].variableIndex(v) << " ";
            }
            std::cout << std::endl;
        }

        if(i < pgm_->Model()->numberOfFactors())
        {
            std::cout << "\tFactor " << i << " has " << (*(pgm_->Model()))[i].numberOfVariables() << " variables in PGM\n\t\t";

            for(size_t v = 0; v < (*(pgm_->Model()))[i].numberOfVariables(); v++)
            {
                std::cout << (*(pgm_->Model()))[i].variableIndex(v) << " ";
            }

            std::cout << std::endl;
        }
    }

    exit(0);
}

void pgmlink::DualDecompositionConservationTracking::configure_hard_constraints(
        const pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::SubGmType& subGM,
        size_t sub_gm_index,
        pgmlink::DualDecompositionConservationTracking::DualDecompositionSubGradient::InfType& optimizer)
{
    assert(hypotheses_graph_ != NULL);

    LOG(logWARNING) << "Found subproblem with " << subGM.numberOfVariables() << " variables";
    current_sub_gm_id_ = sub_gm_index;
    current_sub_optimizer_ = &optimizer;

    // add constraint
    pgmlink::ConservationTracking::add_constraints(*hypotheses_graph_, boost::bind(
                                                       &pgmlink::DualDecompositionConservationTracking::add_constraint,
                                                       this, &optimizer, sub_gm_index, _1, _2, _3, _4, _5, _6));

    // constraint_debug_output(subGM);
}

void pgmlink::DualDecompositionConservationTracking::debug_graph_output(GraphicalModelType* model)
{
    // count factor orders
    std::map<size_t, size_t> factorOrders;

    std::cout << "Model has " << model->numberOfFactors() << " factors and "
              << model->numberOfVariables() << " variables" << std::endl;
    for(size_t i = 0; i < model->numberOfFactors(); i++)
    {
        size_t order = (*model)[i].numberOfVariables();
        if(factorOrders.find(order) == factorOrders.end())
        {
            factorOrders[order] = 1;
        }
        else
            factorOrders[order]++;
    }

    std::cout << "Factor Variable Count:" << std::endl;
    for(std::map<size_t, size_t>::iterator it = factorOrders.begin(); it != factorOrders.end(); ++it)
    {
        std::cout << "\t[" << it->first << "] = " << it->second << std::endl;
    }

    // count variable order
    std::map<size_t, size_t> variableOrders;

    for(size_t i = 0; i < model->numberOfVariables(); i++)
    {
        size_t order = model->numberOfFactors(i);

        if(variableOrders.find(order) == variableOrders.end())
        {
            variableOrders[order] = 1;
        }
        else
            variableOrders[order]++;
    }

    std::cout << "Variable Factor Count:" << std::endl;
    for(std::map<size_t, size_t>::iterator it = variableOrders.begin(); it != variableOrders.end(); ++it)
    {
        std::cout << "\t[" << it->first << "] = " << it->second << std::endl;
    }

    std::cout << "Graph decomposed!" << std::endl;
    exit(0);
}

void pgmlink::DualDecompositionConservationTracking::infer()
{
    if (!with_constraints_) {
        opengm::hdf5::save(optimizer_->graphicalModel(), "./conservationTracking.h5", "conservationTracking");
        throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. "
                                 "The conservation tracking factor graph has been saved to file");
    }

    DualDecompositionSubGradient::Parameter dd_parameter;
    GraphicalModelType* model = pgm_->Model();

    std::cout << "Beginning Graph Decomposition" << std::endl;
    std::cout << "Original Graph had: " << model->numberOfFactors() << " factors and "
              << model->numberOfVariables() << " variables" << std::endl;

    opengm::GraphicalModelDecomposer<GraphicalModelType> decomposer;
    dd_parameter.decomposition_ = decomposer.decomposeIntoClosedBlocks(*model, 2);
    std::cout << "Decomposed into " << dd_parameter.decomposition_.numberOfSubModels() << " submodels" << std::endl;

    for(unsigned int i = 0; i < dd_parameter.decomposition_.numberOfSubModels(); i++)
    {
        std::cout << "\tSubproblem " << i << ": "
                  << dd_parameter.decomposition_.numberOfSubFactors(i) << " factors and "
                  << dd_parameter.decomposition_.numberOfSubVariables(i) << " variables" << std::endl;
    }

    dd_parameter.decomposition_.reorder();
    std::cout << "done reordering, completing... " << std::endl;
    dd_parameter.decomposition_.complete();

    if(!dd_parameter.decomposition_.isValid(*model))
    {
        std::cout << "ERROR: Model decomposition is invalid!!!!!!!!!" << std::endl;
    }
    else
    {
        std::cout << "Model decomposition valid!" << std::endl;
    }

//    debug_graph_output(model);

    dd_parameter.decompositionId_ = DualDecompositionSubGradient::Parameter::MANUAL;
    dd_parameter.maximalDualOrder_ = 4;
    dd_parameter.subPara_.verbose_ = true;
    dd_parameter.subPara_.integerConstraint_ = true;
    dd_parameter.subPara_.epGap_ = ep_gap_;

    dd_optimizer_ = new DualDecompositionSubGradient(*model, dd_parameter,
                        boost::bind(&pgmlink::DualDecompositionConservationTracking::configure_hard_constraints, this, _1, _2, _3));

    DualDecompositionSubGradient::VerboseVisitorType visitor;
    opengm::InferenceTermination status = dd_optimizer_->infer(visitor);

    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated abnormally");
    }
}

void pgmlink::DualDecompositionConservationTracking::extractSolution(
        std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> &solution)
{
    opengm::InferenceTermination status = dd_optimizer_->arg(solution);
    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }

    std::cout << "Dual Decomposition - Found solution:\n";
    for(size_t i = 0; i < solution.size(); i++)
    {
        std::cout << solution[i] << " ";
    }
    std::cout << std::endl;
}

void pgmlink::DualDecompositionConservationTracking::add_constraints(const pgmlink::HypothesesGraph &g)
{
    hypotheses_graph_ = &g;
}

size_t pgmlink::DualDecompositionConservationTracking::cplex_id(size_t opengm_id, size_t state)
{
    // map indices to subproblem indices
    typedef opengm::GraphicalModelDecomposition::SubVariableListType SubVarListType;
    const std::vector<SubVarListType>& sub_variable_list =
            dd_optimizer_->parameter().decomposition_.getVariableLists();

    std::size_t new_id = -1; // init with max id

    for(SubVarListType::const_iterator sub_var_it = sub_variable_list[opengm_id].begin();
        sub_var_it != sub_variable_list[opengm_id].end();
        ++sub_var_it)
    {
        if(sub_var_it->subModelId_ == current_sub_gm_id_)
        {
            // use remapping
            new_id = sub_var_it->subVariableId_;
        }
    }

    if(new_id == (size_t)-1)
    {
        // return an "error" value
        LOG(pgmlink::logDEBUG4) << "OpenGM ID " << opengm_id << " does not exist in subproblem " << current_sub_gm_id_;
        return -1;
    }

    // add constraint if the indices are present in the current subproblem
    return current_sub_optimizer_->lpNodeVi(new_id, state);
}
