#ifndef MULTI_HYPOTHESES_GRAPH
#define MULTI_HYPOTHESES_GRAPH

// stl
#include <string>
#include <vector>

// boost
#include <boost/shared_ptr.hpp>

// lemon
#include <lemon/list_graph.h>

// pgmlink
#include <pgmlink/graph.h>
#include <pgmlink/hypotheses.h>
#include <pgmlink/traxels.h>
#include <pgmlink/region_graph.h>


namespace pgmlink {
/*class TagNode;
  
  class EventNode;
    
  class ObjectNode; */

class Tag;

class ConnectionTag;

class ConflictTag;
  
class MultiHypothesesGraph;

class MultiHypothesesGraphBuilder;

typedef std::vector<boost::shared_ptr<RegionGraph> > RegionGraphVector;

typedef boost::shared_ptr<std::vector<boost::shared_ptr<RegionGraph> > >
RegionGraphVectorPtr;

typedef boost::shared_ptr<MultiHypothesesGraph> MultiHypothesesGraphPtr;

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
  MultiHypothesesGraph();
  unsigned maximum_timestep();
 private:
  unsigned maximum_timestep_;
};


class MultiHypothesesGraphBuilder {
 public:
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
 private:
  void add_nodes(RegionGraphVectorPtr source_graphs,
                 MultiHypothesesGraphPtr dest_graph);
  void add_nodes_at(RegionGraphPtr source_graph,
                    MultiHypothesesGraphPtr dest_graph,
                    unsigned timestep);
  void add_node(RegionGraphPtr source_graph,
                MultiHypothesesGraphPtr dest_graph,
                const RegionGraph::Node& source_node,
                int timestep);
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


    
  
  
  
    
}


#endif /* MULTI_HYPOTHESES_GRAPH */
