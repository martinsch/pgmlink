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


namespace pgmlink {
  class TagNode;
  
  class EventNode;

  class ObjectNode;

  class Tag;

  class ConnectionTag;

  class ConflictTag;
  
  class MultiHypothesesGraph;

  class MultiHypothesesGraphBuilder;

  typedef boost::shared_ptr<std::vector<boost::shared_ptr<RegionGraph> > >
  RegionGraphVector;


  class TagNode : public lemon::ListGraph::Node {
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

  const std::string EventNode::tag = "event";


  class ObjectNode : public TagNode {

  };

  const std::string ObjectNode::tag = "object";
  

  class Tag {
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
  };


  class MultiHypothesesGraph : public PropertyGraph<lemon::ListGraph> {
  public:
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
        max_nearest_neighbors_(max_nearest_neighbors),
        distance_threshold_(distance_threshold),
        forward_backward_(forward_backward),
        consider_divisions_(consider_divisions),
        division_threshold_(division_threshold) {}
      unsigned max_nearest_neighbors_;
      feature_type distance_threshold_;
      bool forward_backward_;
      bool consider_divisions_;
      double division_threshold_;
    };
        
    MultiHypothesesGraphBuilder(const Options& options = Options()) :
      options_(options) {}
    boost::shared_ptr<MultiHypothesesGraph>
    build(RegionGraphVector graph);
  private:
    Options options_;
  };


    
  
  
  
    
}


#endif /* MULTI_HYPOTHESES_GRAPH */
