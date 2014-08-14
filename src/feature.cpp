#include "pgmlink/feature.h"
#include "pgmlink/log.h"
#include <cmath>
#include <string>

using namespace std;

namespace pgmlink {
  ////
  //// class GeometryDivision2
  ////
  double GeometryDivision2::operator()(const Traxel& ancestor,
				       const Traxel& child1,
				       const Traxel& child2) const {
    double feature = 0;    

    if( angle_constraint_ ) {
      double pi = acos(-1.);
      double angle = ancestor.angle(child1, child2);
      if (angle/pi < min_angle_) {
	feature += 10000000000000;
      }
    }
	
    if( distance_dependence_ ) {
      feature += (ancestor.distance_to(child1)-mean_div_dist_) * (ancestor.distance_to(child1)-mean_div_dist_);
      feature += (ancestor.distance_to(child2)-mean_div_dist_) * (ancestor.distance_to(child2)-mean_div_dist_);
    }

    return feature;
  }



  ////
  //// class KasterDivision2
  ////
  double KasterDivision2::operator()(const Traxel& ancestor,
				     const Traxel& child1,
				     const Traxel& child2) const {
    double feature = 0;    
  
    feature += ancestor.distance_to(child1) * ancestor.distance_to(child1);
    feature += ancestor.distance_to(child2) * ancestor.distance_to(child2);
    return w_*(feature + div_cost_);
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

  double get_detection_prob(const Traxel& tr, size_t state) {
	  FeatureMap::const_iterator it = tr.features.find("detProb");	
	  if (it == tr.features.end()) {
		  throw runtime_error("get_detection_prob(): detProb feature not in traxel");
	  }
	  double det_prob = it->second[state];
	  LOG(logDEBUG3) << "get_detection_prob(): " << det_prob;
	  return det_prob;
  }

  double get_division_prob(const Traxel& tr) {
	  FeatureMap::const_iterator it = tr.features.find("divProb");
	  if (it == tr.features.end()) {
		  throw runtime_error("get_division_prob(): divProb feature not in traxel");
	  }
	  double div_prob = it->second[0];
	  LOG(logDEBUG3) << "get_division_prob(): " << div_prob;
	  return div_prob;
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
  
  
////
//// class NegLnDetection
////
double NegLnDetection::operator ()(const Traxel& tr, size_t state) const {
	double arg = get_detection_prob(tr, state);
    if(arg < 0.0000000001) arg = 0.0000000001;
	return w_*-1*log(arg);
}


////
//// class NegLnDivision
////
double NegLnDivision::operator ()(const Traxel& tr, size_t state) const {
	double arg = get_division_prob(tr);
	if (state == 0) {
		arg = 1 - arg;
	}
    if(arg <0.0000000001) arg = 0.0000000001;
	return w_*-1*log(arg);
}


////
//// class NegLnTransition
////
double NegLnTransition::operator ()(const double dist_prob) const {
	double arg = dist_prob;
    if(arg < 0.0000000001) arg = 0.0000000001;
	return w_*-1*log(arg);
}


////
//// class NegLnConstant
////
double NegLnConstant::operator ()(size_t state) const {
	if (state > sizeof(prob_vector_)/sizeof(double)) {
		  throw runtime_error("NegLnConstant(): state must not be larger than the size of the prob. vector");
	}
	double arg = prob_vector_[state];
	LOG(logDEBUG3) << "NegLnConstant(): arg = " << arg;
	if(arg == 0) arg = 0.0000000001;
	return w_*-1*log(arg);
}



  ////
  //// class SquaredDistance
  ////
  double SquaredDistance::operator()(const Traxel& t1,
				     const Traxel& t2) const {
    return t1.distance_to(t2);
  }
    

    
  ////
  //// class KasterDivison
  ////
  double KasterDivision::operator()(const Traxel& ancestor,
				    const Traxel& child1,
				    const Traxel& child2) const {
    SquaredDistance sqd;
    return sqd(ancestor, child1) + sqd(ancestor, child2) + div_cost_;
  }
    


  ////
  //// Class GeometryDivision
  ////
  double GeometryDivision::operator()(const Traxel& ancestor,
				      const Traxel& child1,
				      const Traxel& child2) const {
    double pi = acos(-1.);
    double angle = ancestor.angle(child1, child2);
    double feature = 0;
    if (angle/pi < 0.8) {
      feature = 10000000000000;
    } 
    
    return feature;
  }



  ////
  //// class CellnessDivision
  ////
  CellnessDivision::CellnessDivision(double diffCellness, double absCellness):
    param_diff_c( diffCellness ), param_abs_c( absCellness )
  {}

  double CellnessDivision::operator()(const Traxel& ancestor, 
                                      const Traxel& child1,
                                      const Traxel& child2) const {
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
                                  const Traxel& to) const {
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
  


  double CellnessDisappearance::operator()(const Traxel& from) const {
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



  double CellnessAppearance::operator()(const Traxel& to) const {
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

  ////
  //// SpatialBorderAwareWeight
  ////
 double SpatialBorderAwareWeight::operator()( const Traxel& tr ) const {
	double t = tr.Timestep;
	double x = tr.X(), y = tr.Y(), z = tr.Z();
	double distance_to_border = fov_.spatial_distance_to_border(t,x,y,z,relative_);
	LOG(logDEBUG4) << "SpatialBorderAwareWeight(): distance to border = " << distance_to_border;
	if( distance_to_border < margin_) {
	  double linear_cost = (distance_to_border / margin_) * cost_; //normalize distance within the border to range (0,1)
	  LOG(logDEBUG4) << "SpatialBorderAwareWeight(): distance smaller than margin, energy: " << linear_cost;
	  return linear_cost;
	} else {
		LOG(logDEBUG4) << "SpatialBorderAwareWeight(): distance greater than margin, energy: " << cost_;
		return cost_;
	}
  }


} /* namespace pgmlink */
