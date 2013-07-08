#ifndef MULTI_HYPOTHESES_SEGMENTATION_H
#define MULTI_HYPOTHESES_SEGMENTATION_H

// stl
#include <vector>
#include <set>
#include <map>

// vigra
#include <vigra/multi_array.hxx> // vigra::MultiArray
#include <vigra/multi_iterator.hxx> // vigra::MultiIterator
#include <vigra/multi_iterator_coupled.hxx> // vigra::CoupledIteratorType
#include <vigra/tinyvector.hxx> // vigra::TinyVector

// boost
#include <boost/shared_ptr.hpp>

// lemon
#include <lemon/list_graph.h> // lemon::ListGraph
#include <lemon/maps.h> // lemon::IterableValueMap

// pgmlink
#include <pgmlink/clustering.h>
#include <pgmlink/traxels.h>
#include <pgmlink/region_graph.h>
#include <pgmlink/graph.h>


namespace pgmlink {

  template <typename T, typename U>
  struct AdjacencyList;

  template <typename T, typename U>
  struct AdjacencyListPtr;

  template <typename T, typename U>
  struct PixelActorPtr;

  template <int N, typename T>
  struct NeighborhoodVisitorPtr;

  class MultiSegmenter;
  
  class MultiSegmenterBuilder;

  class MultiSegmentContainer;

  class ConnectedComponentsToMultiSegments;

  template <typename T, typename U>
  class ParentConnectedComponentBase;

  template <int N, typename T>
  class AdjacencyListBuilder;

  template <typename T, typename U>
  class ListInserterBase;
  
  template <int N, typename T>
  class NeighborhoodVisitorBase;

  template <typename T, typename U>
  class PixelActorBase;

  typedef boost::shared_ptr<MultiSegmenter> MultiSegmenterPtr;

  typedef boost::shared_ptr<std::vector<vigra::MultiArray<2, label_type> > > AssignmentListPtr;

  typedef RegionGraph AdjacencyGraph;

  typedef boost::shared_ptr<AdjacencyGraph > AdjacencyGraphPtr;


  ////
  //// AdjacencyList
  ////
  template <typename T, typename U>
  struct AdjacencyList {
    typedef std::map<T, std::set<U> > type;
  };


  ////
  //// AdjacencyListPtr
  ////
  template <typename T, typename U>
  struct AdjacencyListPtr {
    typedef boost::shared_ptr<typename AdjacencyList<T, U>::type > type;
  };


  ////
  //// ListInserterPtr
  ////
  template <typename T, typename U>
  struct ListInserterPtr {
    typedef typename boost::shared_ptr<ListInserterBase<T, U> > type;
  };


  ////
  //// PixelActorPtr
  ////
  template <typename T, typename U>
  struct PixelActorPtr {
    typedef typename boost::shared_ptr<PixelActorBase<T, U> > type;
  };


  ////
  //// NeighborhoodVisitorPtr
  ////
  template <int N, typename T>
  struct NeighborhoodVisitorPtr {
    typedef boost::shared_ptr<NeighborhoodVisitorBase<N, T> > type;
  };


  ////
  //// ParentConnectedComponentPtr
  ////
  template <typename T, typename U>
  struct ParentConnectedComponentPtr {
    typedef boost::shared_ptr<ParentConnectedComponentBase<T, T> > type;
  };


  ////
  //// MultiSegmenter
  ////
  class MultiSegmenter {
  private:
    const std::vector<unsigned>& n_clusters_;
    ClusteringPtr clusterer_;
    MultiSegmenter();
  public:
    MultiSegmenter(const std::vector<unsigned>& n_clusters,
                   ClusteringPtr clusterer);
    vigra::MultiArray<2, unsigned> operator()(uint offset = 0);
    unsigned assign(const arma::vec& sample);
  };


  ////
  //// MultiSegmenterBuilder
  ////
  class MultiSegmenterBuilder {
  private:
    const std::vector<unsigned>& n_clusters_;
    ClusteringBuilderPtr clustering_builder_;
    MultiSegmenterBuilder();
  public:
    MultiSegmenterBuilder(const std::vector<unsigned>& n_clusters,
                          ClusteringBuilderPtr clustering_builder);
    MultiSegmenterPtr build(const feature_array& data);
  };


  ////
  //// MultiSegmentContainer
  ////
  class MultiSegmentContainer {
  private:
    const vigra::MultiArrayView<2, label_type> assignments_;
    const feature_array& coordinates_;
    MultiSegmentContainer();
  public:
    MultiSegmentContainer(vigra::MultiArrayView<2, label_type> assignments,
                          const feature_array& coordinates_);
    int to_images(vigra::MultiArrayView<4, label_type> dest);
  };


  ////
  //// ConnectedComponentsToMultiSegments
  ////
  class ConnectedComponentsToMultiSegments {
  private:
    const std::vector<feature_array >& components_coordinates_;
    const std::vector<unsigned>& n_clusters_;
    ClusteringBuilderPtr clustering_builder_;
    label_type starting_index_;

    AssignmentListPtr assignments_;

    ConnectedComponentsToMultiSegments();
  public:
    ConnectedComponentsToMultiSegments(const std::vector<feature_array >& components_coordinates,
                                       const std::vector<unsigned>& n_clusters,
                                       ClusteringBuilderPtr clustering_builder,
                                       label_type starting_index=1);
    
    AssignmentListPtr get_assignments();
    void to_images(vigra::MultiArray<4, label_type>& dest);
  };


  ////
  //// AdjacencyListBuilder
  ////
  template <int N, typename T>
  class AdjacencyListBuilder {
  public:
    typedef typename vigra::CoupledIteratorType<N, T, T>::type Iterator;
  private:
    vigra::MultiArrayView<N, T> label_image_;
    vigra::MultiArrayView<N, T> connected_component_image_;
    typename NeighborhoodVisitorPtr<N, T>::type neighborhood_accessor_;
    typename ParentConnectedComponentPtr<T, T>::type parent_connected_component_;

    AdjacencyListBuilder();
  public:
    AdjacencyListBuilder(vigra::MultiArrayView<N, T> label_image,
                         vigra::MultiArrayView<N, T> connected_component_image,
                         typename NeighborhoodVisitorPtr<N, T>::type neighborhood_accessor,
                         typename ParentConnectedComponentPtr<T, T>::type parent_connected_component);
    void create_adjacency_list();
  };


  ////
  //// ParentConnectedComponentBase
  ////
  template <typename T, typename U>
  class ParentConnectedComponentBase {
  public:
    virtual ~ParentConnectedComponentBase() {}
    virtual void add_to_connected_component(const T& key,
                                            const U& component) {}
  };


  ////
  //// ParentConnectedComponentGraph
  ////
  template <typename T, typename U>
  class ParentConnectedComponentGraph : public ParentConnectedComponentBase<T, U> {
  private:
    AdjacencyGraphPtr adjacency_graph_;
  public:
    ParentConnectedComponentGraph();
    ParentConnectedComponentGraph(AdjacencyGraphPtr adjacency_graph);
    virtual ~ParentConnectedComponentGraph() {}
    virtual void add_to_connected_component(const T& key,
                                            const U& component);
  };

  ////
  //// ListInserterBase
  ////
  template <typename T, typename U>
  class ListInserterBase {
  private:
  public:
    virtual ~ListInserterBase() {}
    virtual void add_to_list(const T& key,
                             const U& value) = 0;
  };


  ////
  //// ListInserterMap
  ////
  template <typename T, typename U>
  class ListInserterMap : public ListInserterBase<T, U> {
  private:
    typename AdjacencyListPtr<T, U>::type adjacency_list_;
  public:
    ListInserterMap();
    ListInserterMap(typename AdjacencyListPtr<T, U>::type adjacency_list);
    virtual ~ListInserterMap() {}
    virtual void add_to_list(const T& key,
                             const U& value);
  };


  ////
  //// ListInserterGraph
  ////
  template <typename T, typename U>
  class ListInserterGraph : public ListInserterBase<T, U> {
  private:
    AdjacencyGraphPtr adjacency_graph_;
  public:
    ListInserterGraph();
    ListInserterGraph(AdjacencyGraphPtr adjacency_graph);
    virtual ~ListInserterGraph() {};
    virtual void add_to_list(const T& key,
                             const U& value);
  };


  ////
  //// NeighborhoodVisitorBase
  ////
  template <int N, typename T>
  class NeighborhoodVisitorBase {
  private:
  public:
    ~NeighborhoodVisitorBase() {}
    virtual void visit(const vigra::MultiArrayView<N,T> image,
                       const typename AdjacencyListBuilder<N, T>::Iterator pixel) = 0;
  };


  ////
  //// NeighborhoodVisitor2D4Neighborhood
  ////
  template <int N, typename T>
  class NeighborhoodVisitorAdjacentPlusOne : public NeighborhoodVisitorBase<N, T> {
  private:
    typename PixelActorPtr<T, T>::type pixel_actor_;
    NeighborhoodVisitorAdjacentPlusOne();
  public:
    NeighborhoodVisitorAdjacentPlusOne(typename PixelActorPtr<T, T>::type pixel_actor);
    virtual void visit(const vigra::MultiArrayView<N, T> image,
                       const typename AdjacencyListBuilder<N, T>::Iterator pixel);
  };


  ////
  //// PixelActorBase
  ////
  template <typename T, typename U>
  class PixelActorBase {
  private:
  public:
    virtual ~PixelActorBase() {}
    virtual void act(const T& pixel_value,
                     const U& comparison_value) = 0;
  };


  ////
  //// PixelActorFindNeighbors
  ////
  template <typename T>
  class PixelActorFindNeighbors : public PixelActorBase<T, T> {
  private:
    typename ListInserterPtr<T, T>::type inserter_;
    PixelActorFindNeighbors();
  public:
    PixelActorFindNeighbors(typename ListInserterPtr<T, T>::type inserter);
    virtual void act(const T& pixel_value,
                     const T& comparison_value);
  };


  ////
  //// RegionMergingPolicyBase
  ////
  class RegionMergingPolicyBase {
  public:
    virtual ~RegionMergingPolicyBase() {}
    virtual void merge() = 0;
  };


  ////
  //// RegionMergingGraph
  ////
  class RegionMergingGraph : public RegionMergingPolicyBase {
  private:
    AdjacencyGraphPtr adjacency_graph_;
    unsigned maximum_merges_per_connected_component_;
    unsigned maximum_merges_per_patch_;
  public:
    RegionMergingGraph();
    RegionMergingGraph(AdjacencyGraphPtr adjacency_graph,
                       unsigned maximum_merges_per_connected_component,
                       unsigned maximum_merges_per_patch);
    virtual void merge();
  };
  

  /* IMPLEMENTATIONS */


  ////
  //// AdjacencyListBuilder
  ////
  template <int N, typename T>
  AdjacencyListBuilder<N, T>::AdjacencyListBuilder() :
    label_image_(),
    connected_component_image_(),
    neighborhood_accessor_(),
    parent_connected_component_() {

  }


  template <int N, typename T>
  AdjacencyListBuilder<N, T>::AdjacencyListBuilder(vigra::MultiArrayView<N, T> label_image,
                                                   vigra::MultiArrayView<N, T> connected_component_image,
                                                   typename NeighborhoodVisitorPtr<N, T>::type neighborhood_accessor,
                                                   typename ParentConnectedComponentPtr<T, T>::type parent_connected_component) :
    label_image_(label_image),
    connected_component_image_(connected_component_image),
    neighborhood_accessor_(neighborhood_accessor),
    parent_connected_component_(parent_connected_component) {

  }


  template<int N, typename T>
  void AdjacencyListBuilder<N, T>::create_adjacency_list() {
    if (!label_image_.hasData()) {
      return;
    }
    // NeighborhoodVisitorBase<N, label_type>::PixelActorPtr<T, T>::type pixel_actor(new PixelActorFindNeighbors<label_type>);
    Iterator start = createCoupledIterator(label_image_, connected_component_image_);
    Iterator end = start.getEndIterator();
    for (Iterator it = start; it != end; ++it) {
      neighborhood_accessor_->visit(label_image_,
                                    it);
      parent_connected_component_->add_to_connected_component(it.get<1>(), it.get<2>());
    }
  }


  ////
  //// ParentConnectedComponentGraph
  ////
  template <typename T, typename U>
  ParentConnectedComponentGraph<T, U>::ParentConnectedComponentGraph() :
    adjacency_graph_(new AdjacencyGraph) {
    adjacency_graph_->add(node_connected_component());
  }


  template <typename T, typename U>
  ParentConnectedComponentGraph<T, U>::ParentConnectedComponentGraph(AdjacencyGraphPtr adjacency_graph) :
    adjacency_graph_(adjacency_graph) {
    adjacency_graph_->add(node_connected_component());
  }


  template <typename T, typename U>
  void ParentConnectedComponentGraph<T, U>::add_to_connected_component(const T& key,
                                                                       const U& component) {
    AdjacencyGraph::ConnectedComponentMap& connected_component_map =
      adjacency_graph_->get(node_connected_component());
    AdjacencyGraph::LabelMap& label_map = adjacency_graph_->get(node_label());
    AdjacencyGraph::Region region = label_map(key);
    if (region == lemon::INVALID) {
      region = adjacency_graph_->add_region(key);
    }
    connected_component_map.set(region, component);
  }




  ////
  //// NeighborhoodVisitor2DRightLower
  ////
  template <int N, typename T>
  NeighborhoodVisitorAdjacentPlusOne<N, T>::NeighborhoodVisitorAdjacentPlusOne(typename PixelActorPtr<T, T>::type pixel_actor) : 
    pixel_actor_(pixel_actor) {
    
  }



  template <int N, typename T>
  void NeighborhoodVisitorAdjacentPlusOne<N, T>::visit(const vigra::MultiArrayView<N, T> image,
                                                       const typename AdjacencyListBuilder<N, T>::Iterator pixel) {
    vigra::TinyVector<long, N> coordinates = pixel.get<0>();
    for (int dimension = 0; dimension < N; ++dimension) {
      coordinates[dimension] += 1;
      if (image.isInside(coordinates)) {
        pixel_actor_->act(image[coordinates], pixel.get<1>());
      }
      coordinates[dimension] -= 1;
    }
    return;
  }
  

  ////
  //// PixelActorFindNeighbors
  ////
  template <typename T>
  PixelActorFindNeighbors<T>::PixelActorFindNeighbors(typename ListInserterPtr<T, T>::type inserter) :
    inserter_(inserter) {}


  template <typename T>
  void PixelActorFindNeighbors<T>::act(const T& pixel_value,
                                       const T& comparison_value) {
    if (pixel_value != comparison_value &&
        pixel_value != 0 &&
        comparison_value != 0) {
      inserter_->add_to_list(pixel_value, comparison_value);
      inserter_->add_to_list(comparison_value, pixel_value);
    }
  }


  ////
  //// ListInserterMap
  ////
  template <typename T, typename U>
  ListInserterMap<T, U>::ListInserterMap() :
    adjacency_list_(new typename AdjacencyList<T, T>::type) {}

  
  template <typename T, typename U>
  ListInserterMap<T, U>::ListInserterMap(typename AdjacencyListPtr<T, U>::type adjacency_list) :
    adjacency_list_(adjacency_list) {}


  template <typename T, typename U>
  void ListInserterMap<T, U>::add_to_list(const T& key,
                                          const U& value) {
    (*adjacency_list_)[key].insert(value);
  }




  ////
  //// ListInserterGraph
  ////
  template <typename T, typename U>
  ListInserterGraph<T, U>::ListInserterGraph() :
    adjacency_graph_(new AdjacencyGraph) {
    adjacency_graph_->add(node_neighbors())
      .add(arc_dissimilarity())
      .add(node_connected_component());
  }


  template <typename T, typename U>
  ListInserterGraph<T, U>::ListInserterGraph(AdjacencyGraphPtr adjacency_graph) :
    adjacency_graph_(adjacency_graph) {
    adjacency_graph_->add(node_neighbors())
      .add(arc_dissimilarity())
      .add(node_connected_component());
  }


  template <typename T, typename U>
  void ListInserterGraph<T, U>::add_to_list(const T& key,
                                            const U& value) {
    
    AdjacencyGraph::LabelMap& label_map = adjacency_graph_->get(node_label());
    /*AdjacencyGraph::LabelMap::ValueIt value_it = label_map.beginValue();
    for (; value_it != label_map.endValue(); ++value_it) {
      if (*value_it == key) {
        break;
      }
    }*/
    AdjacencyGraph::Region region = label_map(key);
    if (region == lemon::INVALID) {
      region = adjacency_graph_->add_region(key);
    }
    /*if (value_it == label_map.endValue()) {
      region = adjacency_graph_->add_region(key);
    } else {
      region = AdjacencyGraph::LabelMap::ItemIt(label_map, key);
    }
    for (value_it = label_map.beginValue();
         value_it != label_map.endValue(); ++value_it) {
      if (*value_it == value) {
        break;
      }
    } */
    AdjacencyGraph::Region neighbor = label_map(value);
    if (neighbor == lemon::INVALID) {
      neighbor = adjacency_graph_->add_region(value);
    }
    /*if (value_it == label_map.endValue()) {
      neighbor = adjacency_graph_->add_region(value);
    } else {
      neighbor = AdjacencyGraph::LabelMap::ItemIt(label_map, value);
    } */
    adjacency_graph_->get(node_neighbors()).get_value(region).insert(neighbor);
    // also add dissimilarity to arcs?
  }







}

#endif /* MULTI_HYPOTHESES_SEGMENTATION_H */
