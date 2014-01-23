// stl
#include <algorithm>
#include <set>
#include <vector>
#include <cassert>
#include <fstream>

// boost
#include <boost/assert.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>

// lemon
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>

// omp
#include <omp.h>

// pgmlink
#include <pgmlink/multi_hypotheses_graph.h>
#include <pgmlink/traxels.h>
#include <pgmlink/nearest_neighbors.h>
#include <pgmlink/hypotheses.h> // for tags
#include <pgmlink/classifier_auxiliary.h>


namespace pgmlink {
  
  

////
//// MultiHypothesesGraph
////
MultiHypothesesGraph::MultiHypothesesGraph() :
    conflicts_(new ConflictMap) {
  add(node_traxel());


  add(node_timestep());
  add(node_move_features());
  add(node_division_features());
  add(node_count_features());
  add(node_conflict_sets());
  add(arc_from_timestep());
  add(arc_to_timestep());
}


  //
  // write_lgf()
  //
  namespace {
    struct TraxelToStrConverter {
      std::string operator()(const Traxel& t) {
		stringstream ss;
		boost::archive::text_oarchive oa(ss);
		oa & BOOST_SERIALIZATION_NVP(t);
		return ss.str();
      }
    };
    struct VectorTraxelToStrConverter { // vector<Traxel>
	   std::string operator()(const std::vector<Traxel>& t) {
			stringstream ss;
			boost::archive::text_oarchive oa(ss);
			oa & BOOST_SERIALIZATION_NVP(t);
			return ss.str();
		  }
    };
    struct MoveFeaturesToStrConverter { // std::map<Traxel, std::map<Traxel, feature_array> >
    	   std::string operator()(const std::map<Traxel, std::map<Traxel, feature_array> >& t) {
    			stringstream ss;
    			boost::archive::text_oarchive oa(ss);
    			oa & BOOST_SERIALIZATION_NVP(t);
    			return ss.str();
		  }
	};
    struct DivisionFeaturesToStrConverter { // std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> > >
		   std::string operator()(const std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> >& t) {
				stringstream ss;
				boost::archive::text_oarchive oa(ss);
				oa & BOOST_SERIALIZATION_NVP(t);
				return ss.str();
		  }
	};
    struct CountFeaturesToStrConverter { //std::vector<float>
		   std::string operator()(const float& t) {
				stringstream ss;
				boost::archive::text_oarchive oa(ss);
				oa & BOOST_SERIALIZATION_NVP(t);
				return ss.str();
		  }
	};
    struct ConflictSetsToStrConverter { //std::vector<std::vector<unsigned> >
		   std::string operator()(const std::vector<std::vector<unsigned> >& t) {
				stringstream ss;
				boost::archive::text_oarchive oa(ss);
				oa & BOOST_SERIALIZATION_NVP(t);
				return ss.str();
		  }
	};
  }

  void write_lgf( const MultiHypothesesGraph& g, std::ostream& os, bool /*with_n_traxel*/ ) {
	LOG(logDEBUG) << "MultiHypothesesGraph::write_lgf entered";
    lemon::DigraphWriter<HypothesesGraph> writer( g, os );
    writer.
      nodeMap("node_timestep", g.get(node_timestep())).
        // nodeMap("node_regions_in_component", g.get(node_regions_in_component()), VectorTraxelToStrConverter()).
      nodeMap("node_move_features", g.get(node_move_features()), MoveFeaturesToStrConverter()).
      nodeMap("node_division_features", g.get(node_division_features()), DivisionFeaturesToStrConverter()).
      nodeMap("node_count_features", g.get(node_count_features()), CountFeaturesToStrConverter()).
      nodeMap("node_conflict_sets", g.get(node_conflict_sets()), ConflictSetsToStrConverter()).
      arcMap("arc_from_timestep", g.get(arc_from_timestep())).
      arcMap("arc_to_timestep", g.get(arc_to_timestep())).
    //if(with_n_traxel) {
      nodeMap("node_traxel", g.get(node_traxel()), TraxelToStrConverter());
    //}
    writer.run();
  }



//
// read_lgf()
//
namespace {
struct StrToTraxelConverter {
	Traxel operator()(const std::string& s) {
		stringstream ss(s);
		boost::archive::text_iarchive ia(ss);
		Traxel t;
		ia >> t;
		return t;
	}
};

struct StrToVectorTraxelConverter {
	std::vector<Traxel> operator()(const std::string& s) {
		stringstream ss(s);
		boost::archive::text_iarchive ia(ss);
		std::vector<Traxel> t;
		ia >> BOOST_SERIALIZATION_NVP(t);
		return t;
	}
};
struct StrToMoveFeaturesConverter {
	std::map<Traxel, std::map<Traxel, feature_array> > operator()(
			const std::string& s) {
		stringstream ss(s);
		boost::archive::text_iarchive ia(ss);
		std::map<Traxel, std::map<Traxel, feature_array> > t;
		ia >> BOOST_SERIALIZATION_NVP(t);
		return t;
	}
};
struct StrToDivisionFeaturesConverter {
	std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> > operator()(const std::string& s) {
		stringstream ss(s);
		boost::archive::text_iarchive ia(ss);
		std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> > t;
		ia >> BOOST_SERIALIZATION_NVP(t);
		return t;
	}
};
struct StrToCountFeaturesConverter {
	float operator()(const std::string& s) {
		stringstream ss(s);
		boost::archive::text_iarchive ia(ss);
		float t;
		ia >> BOOST_SERIALIZATION_NVP(t);
		return t;
	}
};
struct StrToConflictSetsConverter {
	std::vector<std::vector<unsigned> > operator()(const std::string& s) {
		stringstream ss(s);
		boost::archive::text_iarchive ia(ss);
		std::vector<std::vector<unsigned> >  t;
		ia >> BOOST_SERIALIZATION_NVP(t);
		return t;
	}
};
}

void read_lgf( MultiHypothesesGraph& g, std::istream& is, bool /*with_n_traxel*/ ) {
	LOG(logDEBUG) << "MultiHypothesesGraph::read_lgf entered";
  lemon::DigraphReader<MultiHypothesesGraph> reader( g, is );
  reader.
      nodeMap("node_timestep", g.get(node_timestep())).
      // nodeMap("node_regions_in_component", g.get(node_regions_in_component()),StrToVectorTraxelConverter()).
      nodeMap("node_move_features", g.get(node_move_features()),StrToMoveFeaturesConverter()).
      nodeMap("node_division_features", g.get(node_division_features()), StrToDivisionFeaturesConverter()).
      nodeMap("node_count_features", g.get(node_count_features()), StrToCountFeaturesConverter()).
      nodeMap("node_conflict_sets", g.get(node_conflict_sets()), StrToConflictSetsConverter()).
      arcMap("arc_from_timestep", g.get(arc_from_timestep())).
      arcMap("arc_to_timestep", g.get(arc_to_timestep())).
  	  nodeMap("node_traxel", g.get(node_traxel()), StrToTraxelConverter());
//  if( with_n_traxel ) {
//    g.add(node_traxel());
//    reader.nodeMap("node_traxel", g.get(node_traxel()), StrToTraxelConverter());
//  }
  reader.run();
}

namespace {
void clear_file(const std::string& filename) {
  std::ofstream file;
  file.open(filename.c_str(), std::ios::trunc);
  file.close();
}

}

void MultiHypothesesGraph::add_classifier_features(ClassifierStrategy* move,
                                                   ClassifierStrategy* division,
                                                   ClassifierStrategy* count,
                                                   ClassifierStrategy* detection) {
  // LOG(logDEBUG) << "MultiHypothesesGraph::add_classifier_features: entered";
  // add(node_move_features());
  // add(node_division_features());
  // add(node_count_features());

  // clear_file("classifier_move.log");
  // clear_file("classifier_division.log");
  // clear_file("classifier_detection.log");
  
  // DivisionFeatureMap& division_map = get(node_division_features());
  // MoveFeatureMap& move_map = get(node_move_features());
  // ContainedRegionsMap& regions = get(node_regions_in_component());
  
  // const std::vector<NodeT>& nodes = this->nodes;

  // // allocate space for results before going parallel
  // for (size_t i = 0; i < nodes.size(); ++i) { // need this for omp
  //   	Node n = this->nodeFromId(i);
  //   	assert(this->valid(n));
  //     //  for (NodeIt n(*this); n != lemon::INVALID; ++n) {
  //     LOG(logDEBUG1) << "MultiHypothesesGraph::add_classifier_features: classifying region "
  //                    << get(node_traxel())[n].Id << " at timestep "
  //                    << get(node_traxel())[n].Timestep;
  //     std::vector<Traxel>& sources = regions.get_value(n);
  //     std::vector<Traxel> targets;
  //     for (OutArcIt a(*this, n); a != lemon::INVALID; ++a) {
  //       const std::vector<Traxel>& target = regions[this->target(a)];
  //       LOG(logDEBUG4) << "MultiHypothesesGraph::add_classifier_features: added "
  //                      << target.size() << " regions to target";
  //       targets.insert(targets.end(), target.begin(), target.end());
  //     }
  //     LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: initializing moves";
  //     move->classify(sources, targets, move_map.get_value(n), false); // with_predict = false
  //     LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: initializing divisions";
  //     division->classify(sources, targets, division_map.get_value(n), false); // with_predict = false
  //     LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: initializing detections";
  //     detection->classify(sources, false); // with_predict = false
  //   }


  // // doing the actual predictions in parallel
  // # pragma omp parallel for
  // for (size_t i = 0; i < nodes.size(); ++i) { // need this for omp
  // 	Node n = this->nodeFromId(i);
  // 	assert(this->valid(n));
  //   //  for (NodeIt n(*this); n != lemon::INVALID; ++n) {
  //   LOG(logDEBUG1) << "MultiHypothesesGraph::add_classifier_features: classifying region "
  //                  << get(node_traxel())[n].Id << " at timestep "
  //                  << get(node_traxel())[n].Timestep;
  //   std::vector<Traxel>& sources = regions.get_value(n);
  //   std::vector<Traxel> targets; 
  //   for (OutArcIt a(*this, n); a != lemon::INVALID; ++a) {
  //     const std::vector<Traxel>& target = regions[this->target(a)];
  //     LOG(logDEBUG4) << "MultiHypothesesGraph::add_classifier_features: added "
  //                    << target.size() << " regions to target";
  //     targets.insert(targets.end(), target.begin(), target.end());
  //   }
  //   LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying moves";
  //   move->classify(sources, targets, move_map.get_value(n), true); // with_predict = true
  //   LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying divisions";
  //   division->classify(sources, targets, division_map.get_value(n), true); // with_predict = true
  //   LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying detections";
  //   detection->classify(sources, true); // with_predict = true
  //   LOG(logDEBUG3) << "MultiHypothesesGraph::add_classifier_features: classifying count";
  //   count->classify(sources);
  // }
  
}


void MultiHypothesesGraph::add_conflicts(boost::shared_ptr<std::map<int, std::vector<std::vector<unsigned> > > > conflicts) {
  conflicts_ = conflicts;
}


const ConflictMap& MultiHypothesesGraph::get_conflicts() const {
  return *conflicts_;
}


void MultiHypothesesGraph::remove_traxel_features() {
	// static const string keep_names[] = {"com", "conflicts", "count_prediction", "cardinality", "detProb"};
	// ContainedRegionsMap& regions = get(node_regions_in_component());
	// const std::vector<NodeT>& nodes = this->nodes;

	// for (std::size_t i = 0; i < nodes.size(); ++i) {
	// 	Node n = this->nodeFromId(i);
	// 	assert(this->valid(n));
	// 	LOG(logDEBUG1)
	// 			<< "MultiHypothesesGraph::remove_traxel_features "
	// 			<< get(node_traxel())[n].Id << " at timestep "
	// 			<< get(node_traxel())[n].Timestep;
	// 	std::vector<Traxel>& traxels = regions.get_value(n);

	// 	for(std::vector<Traxel>::iterator trax_it = traxels.begin(); trax_it != traxels.end(); ++trax_it) {
	// 		FeatureMap& feats_old = trax_it->features;
	// 		FeatureMap feats_new;
	// 		for(std::size_t i = 0; i < sizeof(keep_names)/sizeof(keep_names[0]); ++i) {
	// 			FeatureMap::iterator feats_old_it = feats_old.find(keep_names[i]);
	// 			if( feats_old_it != feats_old.end() ) {
	// 				LOG(logDEBUG4) << "keep_names[i]: " << keep_names[i];
	// 				feats_new[keep_names[i]] = feature_array(feats_old_it->second);
	// 			}
	// 		}

	// 		feats_old.clear();

	// 		for(FeatureMap::const_iterator feat_it = feats_new.begin(); feat_it != feats_new.end(); ++feat_it) {
	// 			LOG(logDEBUG4) << feat_it->first << ": " << feat_it->second.size();
	// 			feats_old[feat_it->first] = feat_it->second;
	// 		}
	// 	}
	// }

}


////
//// class SingleTimestepTraxel_MultiHypothesesBuilder
////
HypothesesGraph* SingleTimestepTraxel_MultiHypothesesBuilder::construct() const {
  HypothesesGraph* graph = new MultiHypothesesGraph();
  graph->add(node_traxel());
  return graph;
}


HypothesesGraph* SingleTimestepTraxel_MultiHypothesesBuilder::add_nodes(HypothesesGraph* graph) const {
    LOG(logDEBUG) << "SingleTimestepTraxel_MultiHypothesesBuilder::add_nodes(): entered";

    graph = SingleTimestepTraxel_HypothesesBuilder::add_nodes(graph);

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_m = graph->get(node_traxel());
    MultiHypothesesGraph::ConnectedComponentMap& component_m = graph->get(node_connected_component());

    for(HypothesesGraph::NodeIt node(*graph); node != lemon::INVALID; ++node) {
      const Traxel& trax = traxel_m[node];
      component_m.set(node, std::make_pair(trax.Timestep, trax.Component));
    }
    return graph;
}


} // namespace pgmlink



