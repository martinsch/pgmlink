#ifndef MULTI_HYPOTHESES_GRAPH
#define MULTI_HYPOTHESES_GRAPH

// stl
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <cassert>

// boost
#include <boost/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/export.hpp>

// lemon
#include <lemon/list_graph.h>
#include <lemon/maps.h>

// pgmlink
#include "pgmlink/graph.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"
#include "pgmlink/region_graph.h"
// #include "pgmlink/classifier_auxiliary.h"



namespace pgmlink {
/*class TagNode;
  
  class EventNode;
    
  class ObjectNode; */

class Tag;

class ConnectionTag;

class ConflictTag;
  
class MultiHypothesesGraph;

class MultiHypothesesGraphBuilder;

class MultiHypothesesTraxelStore;

// class MultiHypothesesTraxelStoreBuilder;

class ClassifierStrategy;

typedef std::vector<boost::shared_ptr<RegionGraph> > RegionGraphVector;

typedef boost::shared_ptr<std::vector<boost::shared_ptr<RegionGraph> > >
RegionGraphVectorPtr;

typedef boost::shared_ptr<MultiHypothesesGraph> MultiHypothesesGraphPtr;

typedef std::map<int, std::map<unsigned, std::vector<Traxel> > > TimestepRegionMap;

typedef std::map<unsigned, std::vector<std::vector<unsigned> > > ConflictSetMap;

typedef std::map<int, std::vector<std::vector<int> > > SignedConflictSetMap;

typedef std::map<int, ConflictSetMap> TimestepConflictSetMap;

template <typename PropertyTag, typename Graph>
struct PropertyValue {
  typedef typename property_map<PropertyTag, Graph>::type::Value type;
};

template <typename PropertyTag, typename Graph>
struct PropertyValueVector {
  typedef std::vector<typename PropertyValue<PropertyTag, Graph>::type> type;
};

template <typename PropertyTag, typename Graph>
struct PropertyValueVectorPtr {
  typedef boost::shared_ptr<typename PropertyValueVector<PropertyTag, Graph>::type> type;
};

////
//// node_regions_in_component
////
struct node_regions_in_component{};

template <typename Graph>
struct property_map<node_regions_in_component, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::vector<Traxel> > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_regions_in_component, Graph>::name = "node_regions_in_component";


////
//// node_move_features
////
struct node_move_features{};

template <typename Graph>
struct property_map<node_move_features, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::map<Traxel, std::map<Traxel, feature_array> > > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_move_features, Graph>::name = "node_move_features";


////
//// node_division_features
////
struct node_division_features{};

template <typename Graph>
struct property_map<node_division_features, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::map<Traxel, std::map<std::pair<Traxel, Traxel>, feature_array> > > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_division_features, Graph>::name = "node_division_features";


////
//// node_count_features
////
struct node_count_features{};

template <typename Graph>
struct property_map<node_count_features, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, feature_type> type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_count_features, Graph>::name = "node_count_features";


////
//// node_conflict_sets
////
struct node_conflict_sets{};

template <typename Graph>
struct property_map<node_conflict_sets, Graph> {
  typedef IterableEditableValueMap<Graph, typename Graph::Node, std::vector<std::vector<unsigned> > > type;
  static const std::string name;
};

template <typename Graph>
const std::string property_map<node_conflict_sets, Graph>::name = "node_conflict_sets";



/* class TagNode : public lemon::ListGraph::Node {
   public:
   static const std::string tag;
   };

   const std::string TagNode::tag = "";


   class EventNode : public TagNode {
   public:
   enum EventType {Move, Division, Appearance, Disappearance, Void};
    
   EventNode(EventType type, unsigned from, unsigned to);

   EventType get_type();
   unsigned from();
   unsigned to();
   private:
   EventType type_;
   unsigned from_;
   unsigned to_;
   feature_array distances_;
   };


   class ObjectNode : public TagNode {

   };

   const std::string ObjectNode::tag = "object"; */
  

/* class Tag {
   public:
   virtual std::string tag() {return ""} const;
   virtual ~Tag();
   };

  
   class ConnectionTag : public Tag {
   public:
   virtual std::string tag() const {return "connection"};
   virtual ~ConnectionTag();
   };

  
   class ConflictTag : public Tag {
   public:
   virtual std::string tag() const {return "conflict"};
   virtual ~ConflictTag();
   };
  
  
   class TaggedArc : public lemon::ListGraph::Arc {
   public:
   TaggedArc(Tag* tag) :
   tag(tag->tag()) {}
   const std::string tag;
   }; */


class MultiHypothesesGraph : public HypothesesGraph {
 public:
  enum EventType {Object, Move, Division, Appearance, Disappearance};
  enum ArcType {Connection, Conflict};
  typedef property_map<node_regions_in_component, base_graph>::type ContainedRegionsMap;
  typedef property_map<node_traxel, base_graph>::type TraxelMap;
  typedef property_map<node_division_features, base_graph>::type DivisionFeatureMap;
  typedef property_map<node_move_features, base_graph>::type MoveFeatureMap;
  typedef property_map<node_count_features, base_graph>::type CountFeatureMap;
  typedef property_map<node_conflict_sets, base_graph>::type ConflictSetMap;


  MultiHypothesesGraph();

  template <typename PropertyTag>
  typename PropertyValueVectorPtr<PropertyTag, base_graph>::type
  get_properties_at(int timestep);

  void add_classifier_features(ClassifierStrategy* move,
                               ClassifierStrategy* division,
                               ClassifierStrategy* count,
                               ClassifierStrategy* detection);

  
 private:
  // boost serialize
  friend class boost::serialization::access;
  template< typename Archive >
	void save( Archive&, const unsigned int /*version*/ ) const;
  template< typename Archive >
	void load( Archive&, const unsigned int /*version*/ );
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  unsigned maximum_timestep_;
};

// lemon graph format (lgf) serialization
  PGMLINK_EXPORT void write_lgf( const MultiHypothesesGraph&, std::ostream& os=std::cout,
		  bool with_n_traxel=false );
  PGMLINK_EXPORT void read_lgf( MultiHypothesesGraph&, std::istream& is=std::cin,
		 bool with_n_traxel=false);


  template< typename Archive >
    void MultiHypothesesGraph::save( Archive& ar, const unsigned int /*version*/ ) const {
	LOG(logDEBUG) << "MultiHypothesesGraph::save entered";
    ar & maximum_timestep_;
    ar & timesteps_;

    bool with_n_traxel = false;
    try {
      this->get(node_traxel());
      with_n_traxel = true;
    } catch( std::runtime_error& e ) {}

    ar & with_n_traxel;

    std::string lgf;
    {
      std::stringstream ss;
      write_lgf(*this, ss, with_n_traxel);
      lgf = ss.str();
    }
    ar & lgf;
  }

  template< typename Archive >
    void MultiHypothesesGraph::load( Archive& ar, const unsigned int /*version*/ ) {
	  LOG(logDEBUG) << "MultiHypothesesGraph::load entered";
    ar & maximum_timestep_;
    ar & timesteps_;

    bool with_n_traxel;
    ar & with_n_traxel;

    std::string lgf;
    ar & lgf;
    {
      std::stringstream ss(lgf);
      read_lgf(*this, ss, with_n_traxel);
    }
   }



class MultiHypothesesGraphBuilder {
 public:
  enum DIRECTION {FORWARD = 1, BACKWARD = -1};
  struct Options {
    Options(unsigned max_nearest_neighbors=2,
            feature_type distance_threshold=50,
            bool forward_backward=false,
            bool consider_divisions=true,
            double division_threshold=0.5) :
        max_nearest_neighbors(max_nearest_neighbors),
        distance_threshold(distance_threshold),
        forward_backward(forward_backward),
        consider_divisions(consider_divisions),
        division_threshold(division_threshold) {}
    unsigned max_nearest_neighbors;
    feature_type distance_threshold;
    bool forward_backward;
    bool consider_divisions;
    double division_threshold;
  };
        
  MultiHypothesesGraphBuilder(const Options& options = Options());
  MultiHypothesesGraphPtr build(RegionGraphVectorPtr graphs);
  MultiHypothesesGraphPtr build(const MultiHypothesesTraxelStore& ts);
 private:
  void add_nodes(RegionGraphVectorPtr source_graphs,
                 MultiHypothesesGraphPtr dest_graph);
  
  void add_nodes(const MultiHypothesesTraxelStore& ts,
                 MultiHypothesesGraphPtr dest_graph);
  
  void add_nodes_at(RegionGraphPtr source_graph,
                    MultiHypothesesGraphPtr dest_graph,
                    unsigned timestep);

  void add_node(RegionGraphPtr source_graph,
                MultiHypothesesGraphPtr dest_graph,
                const RegionGraph::Node& source_node,
                int timestep);

 

  void add_edges(MultiHypothesesGraphPtr graph);
  void add_edges_at(MultiHypothesesGraphPtr graph,
                    int timestep,
                    DIRECTION direction);
  void add_edges_for_node(MultiHypothesesGraphPtr graph,
                          const MultiHypothesesGraph::Node& node,
                          std::map<unsigned, double>& neighbors,
                          int timestep,
                          DIRECTION direction);
  TraxelVectorPtr extract_traxels(RegionGraphPtr graph,
                                  unsigned cc_label);
  TraxelVectorPtr reduce_to_nearest_neighbors(TraxelVectorPtr traxels,
                                              std::map<unsigned, double>& neighbors);
  void create_events_for_component(const Traxel& trax,
                                   TraxelVectorPtr traxels,
                                   RegionGraphPtr single_timestep_graph,
                                   MultiHypothesesGraphPtr graph);
  void add_move_events(const Traxel& trax,
                       std::map<unsigned, double>& neighbors,
                       MultiHypothesesGraphPtr graph);
  void add_division_events(const Traxel& trax,
                           std::map<unsigned, double>& neighbors,
                           MultiHypothesesGraphPtr graph);
  void add_appearance_events(const Traxel& trax,
                             MultiHypothesesGraphPtr graph);
  void add_disappearance_events(const Traxel& trax,
                                MultiHypothesesGraphPtr graph);

  void transfer_nodes(RegionGraphPtr region_graph,
                      MultiHypothesesGraphPtr multi_hypotheses_graph);
    
  Options options_;
  std::map<RegionGraph::Node, MultiHypothesesGraph::Node> reference_map_;
  std::map<MultiHypothesesGraph::Node, RegionGraph::Node> cross_reference_map_;
};


////
//// class MultiHypothesesTraxelStore
////
struct MultiHypothesesTraxelStore {
 public:
  MultiHypothesesTraxelStore() {}

  MultiHypothesesTraxelStore(TimestepRegionMap map, TimestepConflictSetMap conflicts_by_timestep):
	  map(map), conflicts_by_timestep(conflicts_by_timestep) {}

  void add(const Traxel& trax, unsigned component_id);
  void add_conflict_map(int timestep, const ConflictSetMap& conflicts);
  void add_signed_conflict_map(int timestep, const SignedConflictSetMap& conflicts);
  void start_component(const Traxel& trax);
  void make_sane();
  std::string print();

  TimestepRegionMap map;
  TimestepConflictSetMap conflicts_by_timestep;

 private:
  friend class boost::serialization::access;
  template< typename Archive >
   void serialize( Archive&, const unsigned int /*version*/ );
};


template< typename Archive >
void MultiHypothesesTraxelStore::serialize( Archive& ar, const unsigned int /*version*/ ) {
   ar & map;
   ar & conflicts_by_timestep;
}


////
//// class MultiHypothesesTraxelStoreBuilder
////

/* class MultiHypothesesTraxelStoreBuilder {
 public:
  template <int N, typename LABEL_TYPE>
  void build(MultiHypothesesTraxelStore& ts,
             vigra::MultiArrayView<N+1, LABEL_TYPE> arr,
             vigra::MultiArrayView<N, LABEL_TYPE> components,
             unsigned object_layer,
             const Traxel& trax
             );
 private:
  template <int N, typename LABEL_TYPE>
  Traxel& assign_component(MultiHypothesesTraxelStore& ts,
                           vigra::MultiArrayView<N, LABEL_TYPE> arr,
                           vigra::MultiArrayView<N, LABEL_TYPE> components,
                           const Traxel& trax
                           );
}; */

 
/* IMPLEMENTATIONS */



template <typename PropertyTag>
typename PropertyValueVectorPtr<PropertyTag, MultiHypothesesGraph::base_graph>::type
MultiHypothesesGraph::get_properties_at(int timestep) {
  typename PropertyValueVectorPtr<PropertyTag, base_graph>::type
      properties(new typename PropertyValueVector<PropertyTag, base_graph>::type);
  typename property_map<PropertyTag, base_graph>::type& map = get(PropertyTag());
  node_timestep_map& time_map = get(node_timestep());
  for (node_timestep_map::ItemIt it(time_map, timestep);
       it != lemon::INVALID;
       ++it) {
    properties->push_back(map[it]);
  }
  return properties;
}


////
//// class MultiHypothesesTraxelStoreBuilder
////
/* template <int N, typename LABEL_TYPE>
void MultiHypothesesTraxelStoreBuilder::build(MultiHypothesesTraxelStore& ts,
                                              vigra::MultiArrayView<N+1, LABEL_TYPE> arr,
                                              vigra::MultiArrayView<N, LABEL_TYPE> components,
                                              unsigned object_layer,
                                              const Traxel& trax
                                              ) {
  typedef typename vigra::CoupledIteratorType<N>::type ITERATOR;
  unsigned bind_axis = N;
  LABEL_TYPE object_label = trax.Id;

  LOG(logDEBUG4) << "MultiHypothesesTraxelStoreBuilder::build - arr.shape(): " << arr.shape()
                 << ", components.shape(): " << components.shape();
  unsigned number_of_layers = arr.shape()[arr.shape().size()-1];
  assert(object_layer < number_of_layers);
  Traxel& traxel = assign_component<N, LABEL_TYPE>(ts, arr.bindAt(bind_axis, object_layer), components, trax);
  feature_array& conflicts = traxel.features["conflicts"];

  ITERATOR start = createCoupledIterator(arr.bindAt(bind_axis, object_label).shape());
  ITERATOR end = start.getEndIterator();
  for (; start != end; ++start) {
    const vigra::TinyVector<vigra::MultiArrayIndex, N> coords = start.get<0>();
    if (arr.bindAt(bind_axis, object_layer)[coords] == object_label) {
      for (unsigned layer = 0; layer < number_of_layers; ++layer) {
        if (layer != object_layer) {
          LABEL_TYPE label = arr.bindAt(bind_axis, layer)[coords];
          if (label != 0 &&
              std::find(conflicts.begin(), conflicts.end(), label) == conflicts.end()) {
            conflicts.push_back(label);
          }
        }
      }
    }
  }
}


template <int N, typename LABEL_TYPE>
Traxel& MultiHypothesesTraxelStoreBuilder::assign_component(MultiHypothesesTraxelStore& ts,
                                                            vigra::MultiArrayView<N, LABEL_TYPE> arr,
                                                            vigra::MultiArrayView<N, LABEL_TYPE> connected_components,
                                                            const Traxel& trax
                                                            ) {
  LOG(logDEBUG4) << "MultiHypothesesTraxelStoreBuilder::assign_component - arr: " << arr.shape()
                 << ", com: " << connected_components.shape();
  typedef typename vigra::CoupledIteratorType<N, LABEL_TYPE, LABEL_TYPE>::type ITERATOR;
  ITERATOR start = createCoupledIterator(arr, connected_components);
  ITERATOR end = start.getEndIterator();
  int timestep = trax.Timestep;
  LABEL_TYPE object_label = trax.Id;

  LABEL_TYPE component_label = 0;
  std::map<unsigned, std::pair<Traxel, std::vector<Traxel> > >& components = ts.map[timestep];
  for (; start != end; ++start) {
    if (start.get<1>() == object_label) {
      component_label = start.get<2>();
      if (components.find(component_label) == components.end()) {
        ts.start_component(Traxel(component_label, timestep));
      }
      components[component_label].second.push_back(trax);
      components[component_label].second.begin()->features["conflicts"].push_back(object_label);
      components[component_label].second.rbegin()->features["conflicts"].push_back(component_label);
      break;
    }
  }
  return *(components[component_label].second.rbegin());
} */





  
    
}


#endif /* MULTI_HYPOTHESES_GRAPH */




