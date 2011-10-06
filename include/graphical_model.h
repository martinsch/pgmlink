#ifndef GRAPHICAL_MODEL_H
#define GRAPHICAL_MODEL_H

#include <vector>
#include <utility>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <opengm/adder.hxx>
#include <opengm/discretespace.hxx>
#include <opengm/graphicalmodel.hxx>
#include <opengm/explicitfactor.hxx>
#include <opengm/inference/lpcplex.hxx>
#include <vigra/multi_array.hxx>

//#include "hypotheses.h"
#include "graph.h"

using boost::vecS;
using boost::bidirectionalS;
using boost::shared_ptr;

namespace Tracking {

class OpengmMrf {
    public:
    typedef double Energy;
    typedef opengm::DiscreteSpace ogmSpace;
    typedef opengm::ExplicitFactor<Energy> ogmFactor;
    typedef opengm::GraphicalModel<ogmFactor, opengm::Adder> ogmGraphicalModel;
    typedef opengm::Minimizer ogmAccumulator;
    typedef opengm::Inference<ogmGraphicalModel, ogmAccumulator> ogmInference;
    
    OpengmMrf();
    ~OpengmMrf();

    ogmSpace* Space() {return space_; }
    ogmGraphicalModel* Model() {return model_; }
 
    private:
    OpengmMrf(const OpengmMrf&);
    OpengmMrf& operator=(const OpengmMrf&);

    ogmSpace* space_;
    ogmGraphicalModel* model_;
};

} /* namespace Tracking */


#endif /* GRAPHICAL_MODEL_H */
