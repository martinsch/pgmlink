#include "energy.h"
#include "log.h"
#include <cmath>
#include <string>

using namespace std;

namespace Tracking {
////
//// class DivisionEner
////
double GeometryDivision2::operator()(const Traxel& ancestor,
	const Traxel& child1,
	const Traxel& child2) const {
	double energy = 0;    

	if( angle_constraint_ ) {
	  double pi = acos(-1);
	  double angle = ancestor.angle(child1, child2);
	  if (angle/pi < min_angle_) {
	    energy += 10000000000000;
	  }
	}
	
	if( distance_dependence_ ) {
	  energy += (ancestor.distance_to(child1)-mean_div_dist_) * (ancestor.distance_to(child1)-mean_div_dist_);
	  energy += (ancestor.distance_to(child2)-mean_div_dist_) * (ancestor.distance_to(child2)-mean_div_dist_);
	}

	return energy;
}



////
//// class KasterDivision2
////
double KasterDivision2::operator()(const Traxel& ancestor,
	const Traxel& child1,
	const Traxel& child2) const {
  double energy = 0;    
  
  energy += ancestor.distance_to(child1) * ancestor.distance_to(child1);
  energy += ancestor.distance_to(child2) * ancestor.distance_to(child2);
  return w_*(energy + div_cost_);
}

namespace {
  double get_cellness(const Traxel& tr) {
        FeatureMap::const_iterator it = tr.features.find("cellness");
        if(it == tr.features.end()) {
	  throw runtime_error("get_cellness(): cellness feature not in traxel");
	}
        double cellness = it->second[0];
	LOG(logDEBUG3) << "get_cellness(): " << cellness;
	return cellness;
  }
}


////
//// class NegLnCellness
////
double NegLnCellness::operator()(const Traxel& tr) const {
  double cellness = get_cellness(tr);
  if(cellness == 0) cellness = 0.00001;
  return w_*-1*log(cellness);
}


////
//// class NegLnOneMinusCellness
////
double NegLnOneMinusCellness::operator()(const Traxel& tr) const {
  double arg = 1 - get_cellness(tr);
  if(arg == 0) arg = 0.00001;
  return w_*-1*log(arg);
}



    NullaryEnergy::~NullaryEnergy() {}
    UnaryEnergy::~UnaryEnergy() {}
    BinaryEnergy::~BinaryEnergy() {}
    TertiaryEnergy::~TertiaryEnergy() {}

    ConstantEnergy::ConstantEnergy( double energy ) : theEnergy( energy ) {} 
  double ConstantEnergy::operator()(const Traxel& /*t1*/,
				    const Traxel& /*t2*/,
				    const Traxel& /*t3*/,
				    const Traxels& /*prev*/,
				    const Traxels& /*curr*/) const {
	return theEnergy;
    }

  double ConstantEnergy::operator()(const Traxel& /*t1*/,
				    const Traxel& /*t2*/,
				    const Traxels& /*prev*/,
				    const Traxels& /*curr*/) const {
	return theEnergy;
    }

  double ConstantEnergy::operator()(const Traxel& /*t*/,
				    const Traxels& /*prev*/,
				    const Traxels& /*curr*/) const {
	return theEnergy;
    }


  double ConstantEnergy::operator()(const Traxels& /*prev*/,
				    const Traxels& /*curr*/) const {
	 return theEnergy;
    }

    

    ////
    //// class SquaredDistance
    ////
    double SquaredDistance::operator()(const Traxel& t1,
	const Traxel& t2,
				       const Traxels& /*prev*/,
				       const Traxels& /*curr*/) const {
	double dx = abs(t1.features.find(loc_feat_)->second[0] - t2.features.find(loc_feat_)->second[0]);
	double dy = abs(t1.features.find(loc_feat_)->second[1] - t2.features.find(loc_feat_)->second[1]);
	double dz = abs(t1.features.find(loc_feat_)->second[2] - t2.features.find(loc_feat_)->second[2]);
    
	return dx*dx + dy*dy + dz*dz;
    }
    

    
    ////
    //// class KasterDivison
    ////
    double KasterDivision::operator()(const Traxel& ancestor,
                          const Traxel& child1,
                          const Traxel& child2,
				      const Traxels& prev,
				      const Traxels& curr) const {
	SquaredDistance sqd;
	return sqd(ancestor, child1, prev, curr) + sqd(ancestor, child2, prev, curr) + div_cost_;
    }
    


    ////
    //// Class GeometryDivision
    ////
    double GeometryDivision::operator()(const Traxel& ancestor,
                          const Traxel& child1,
                          const Traxel& child2,
					const Traxels& /*prev*/,
					const Traxels& /*curr*/) const {
	double pi = acos(-1);
	double angle = ancestor.angle(child1, child2);
	double energy = 0;
	if (angle/pi < 0.8) {
	    energy = 10000000000000;
	} 

	return energy;
    }



    ////
    //// class CellnessDivision
    ////
    CellnessDivision::CellnessDivision(double diffCellness, double absCellness):
            param_diff_c( diffCellness ), param_abs_c( absCellness )
    {}



    double CellnessDivision::operator()(const Traxel& ancestor, 
                                      const Traxel& child1,
                                      const Traxel& child2,
					const Traxels& /*prev*/,
					const Traxels& /*curr*/) const
    {
        // calculate distance
        FeatureMap::const_iterator pAncestor = ancestor.features.find("com");
        FeatureMap::const_iterator pChild1 = child1.features.find("com");
        FeatureMap::const_iterator pChild2 = child2.features.find("com");

        double dist1 = 0;
        double dist2 = 0;
        if(pAncestor != ancestor.features.end() &&
           pChild1 != child1.features.end() &&
           pChild2 != child2.features.end() )
        {
            //squared distance
            dist1 = std::pow((*pAncestor).second[0]-(*pChild1).second[0],2) +
                    std::pow((*pAncestor).second[1]-(*pChild1).second[1],2) +
                    std::pow((*pAncestor).second[2]-(*pChild1).second[2],2);

            dist2 = std::pow((*pAncestor).second[0]-(*pChild2).second[0],2) +
                    std::pow((*pAncestor).second[1]-(*pChild2).second[1],2) +
                    std::pow((*pAncestor).second[2]-(*pChild2).second[2],2);
        }

        // incorporate cellness
        FeatureMap::const_iterator cAncestor = ancestor.features.find("cellness");
        FeatureMap::const_iterator cChild1 = child1.features.find("cellness");
        FeatureMap::const_iterator cChild2 = child2.features.find("cellness");

        double cellness = 0; // cellness of ancestor
        double dcellness12 = 0; // cellness difference child1-child2

        //std::cout << "From " << ancestor.ID << " to " << child1.ID << " and " << child2.ID << "\n" ;

        if(cAncestor != ancestor.features.end() &&
           cChild1 != child1.features.end() &&
           cChild2 != child2.features.end() )
        {
            // absolute value for ancestor cell
            cellness = (*cAncestor).second[0];

            // squared difference
            dcellness12 = std::pow((*cChild1).second[0]-(*cChild2).second[0],2);

            //std::cout << "Cellness parent " << (*cAncestor).second[0] << " \n" ;
            //std::cout << "Cellness children " << (*cChild1).second[0] << " , " << (*cChild2).second[0] << "\n" ;
        }


        // calculate costs of each term

        // define an epsilon to prevent error at cellness = 0
        double e = 0.0000001;

        // scale the distance
        double cost_dist1 = dist1;
        double cost_dist2 = dist2;

        // penalize a low cellness
        double cost_cellness = - param_abs_c *std::log((cellness*cellness+e)/(1.+e));

        // penalize a big cellness difference of daughter cells
        double cost_dcellness12 = - param_diff_c *std::log((1-dcellness12+e)/(1.+e));

        return cost_dist1 + cost_dist2 + cost_cellness + cost_dcellness12;
    }




    ////
    //// class CellnessMove
    ////
    CellnessMove::CellnessMove(double diffCellness):
            param_diff_c( diffCellness )
    {}



    double CellnessMove::operator()(const Traxel& from, 
                                  const Traxel& to,
				    const Traxels& /*prev*/,
				    const Traxels& /*curr*/) const
    {
        // calculate distance
        FeatureMap::const_iterator pFrom = from.features.find("com");
        FeatureMap::const_iterator pTo = to.features.find("com");

        double dist = 0;
        if(pFrom != from.features.end() && pTo != to.features.end())
        {
            //squared distance
            dist = std::pow((*pFrom).second[0]-(*pTo).second[0],2) +
                   std::pow((*pFrom).second[1]-(*pTo).second[1],2) +
                   std::pow((*pFrom).second[2]-(*pTo).second[2],2);
        }

        // incorporate cellness
        FeatureMap::const_iterator cFrom = from.features.find("cellness");
        FeatureMap::const_iterator cTo = to.features.find("cellness");

        double dcellness = 0;
        if(cFrom != from.features.end() && cTo != to.features.end())
        {
            //squared difference
            dcellness = std::pow((*cFrom).second[0]-(*cTo).second[0],2);
        }

        // calculate costs for each term

        // take distance "as is"
        double cost_dist = dist;

        // define an epsilon to prevent error at log(0)
        double e = 0.0000001;

        //adjust cellness term accordingly
        double cost_dcellness = - param_diff_c *std::log((1-dcellness+e)/(1.+e));

        //std::cout << "From " << from.ID << " to " << to.ID << "\n" ;

        return cost_dist + cost_dcellness;
    }




    CellnessDisappearance::CellnessDisappearance(double disappWeight):
            param_abs_c(disappWeight)
    {}



    double CellnessDisappearance::operator()(const Traxel& from,
                                             const Traxels& /*prev*/,
                                             const Traxels& /*curr*/) const
    {
        // incorporate cellness
        FeatureMap::const_iterator cFrom = from.features.find("cellness");
        if(cFrom == from.features.end()) {
	  throw runtime_error("CellnessDisappearance::operator(): cellness feature not in traxel");
	}

        // define an epsilon to prevent error at cellness = 0
        double e = 0.0000001;

        double cellness = 0;
        if(cFrom != from.features.end())
        {
            //absolute cellness
            cellness = (*cFrom).second[0];
        }
	LOG(logDEBUG4) << "CellnessDisappearance::operator(): cellness " << cellness;
        double cost_cellness = - param_abs_c *std::log((1-cellness+e)/(1.+e));
	LOG(logDEBUG4) << "CellnessDisappearance::operator(): cellness cost " << cost_cellness;
        return cost_cellness;
    }




    CellnessAppearance::CellnessAppearance(double appWeight):
            param_abs_c(appWeight)
    {}



    double CellnessAppearance::operator()(const Traxel& to, 
					  const Traxels& /*prev*/,
					  const Traxels& /*curr*/) const
    {
        // incorporate cellness
        FeatureMap::const_iterator cTo = to.features.find("cellness");
        if(cTo == to.features.end()) {
	  throw runtime_error("CellnessDisappearance::operator(): cellness feature not in traxel");
	}

        // define an epsilon to prevent error at cellness = 0
        double e = 0.00001;

        double cellness = 0;
        if(cTo != to.features.end())
        {
            //absolute cellness
            cellness = (*cTo).second[0];
        }
	LOG(logDEBUG4) << "CellnessAppearance::operator(): cellness " << cellness;
        double cost_cellness = - param_abs_c *std::log((1-cellness+e)/(1.+e));
	LOG(logDEBUG4) << "CellnessAppearance::operator(): cellness cost " << cost_cellness;

        return cost_cellness;
    }


} /* namespace Tracking */
