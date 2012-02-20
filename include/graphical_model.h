#ifndef GRAPHICAL_MODEL_H
#define GRAPHICAL_MODEL_H

#include <vector>
#include <utility>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <opengm/operations/adder.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/lpcplex.hxx>
#include <vigra/multi_array.hxx>

#include "hypotheses.h"
#include "graph.h"

using boost::vecS;
using boost::bidirectionalS;
using boost::shared_ptr;

namespace Tracking {

class OpengmMrf {
   public:
   typedef double Energy;
   typedef opengm::GraphicalModel<Energy, opengm::Adder> ogmGraphicalModel;
   typedef opengm::Factor<ogmGraphicalModel> ogmFactor;
   typedef opengm::Minimizer ogmAccumulator;
   typedef opengm::Inference<ogmGraphicalModel, ogmAccumulator> ogmInference;
   typedef opengm::meta::TypeAtTypeList<ogmGraphicalModel::FunctionTypeList, 0>::type ExplicitFunctionType;
   typedef ogmGraphicalModel::FunctionIdentifier FunctionIdentifier;
    
   OpengmMrf();
   ~OpengmMrf();

   ogmGraphicalModel* Model() {return model_; }
   const ogmGraphicalModel* Model() const {return model_; }
 
   private:
   OpengmMrf(const OpengmMrf&);
   OpengmMrf& operator=(const OpengmMrf&);
   
   ogmGraphicalModel* model_;
};

} /* namespace Tracking */

#endif /* GRAPHICAL_MODEL_H */
