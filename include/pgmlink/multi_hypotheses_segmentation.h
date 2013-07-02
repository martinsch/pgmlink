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

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include <pgmlink/clustering.h>
#include <pgmlink/traxels.h>
#include <pgmlink/graph.h>

// lemon
#include <lemon/list_graph.h> // lemon::ListGraph
#include <lemon/maps.h> // lemon::IterableValueMap


namespace pgmlink {

  template <typename T, typename U>
  struct AdjacencyList;

  template <typename T, typename U>
  struct AdjacencyListPtr;

  class MultiSegmenter;
  
  class MultiSegmenterBuilder;

  class MultiSegmentContainer;

  class ConnectedComponentsToMultiSegments;

  template <int N>
  class AdjacencyListBuilder;

  template <typename T, typename U>
  class ListInserterBase;
  
  template <int N, typename T>
  class NeighborhoodVisitorBase;

  template <typename T, typename U, typename V>
  class PixelActorBase;

  typedef unsigned label_type;

  typedef boost::shared_ptr<MultiSegmenter> MultiSegmenterPtr;

  typedef boost::shared_ptr<std::vector<vigra::MultiArray<2, label_type> > > AssignmentListPtr;

  typedef PropertyGraph<lemon::ListGraph> AdjacencyGraph;

  typedef boost::shared_ptr<AdjacencyGraph > AdjacencyGraphPtr;

  
  ////
  //// node_neighbors
  ////  
  struct node_neighbors {};
  template <typename Graph>
  struct property_map<node_neighbors, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, AdjacencyList<label_type, label_type> > type;
    static const std::string name;
  };


  ////
  //// arc_dissimilarity
  ////
  struct arc_dissimilarity {};
  template <typename Graph>
  struct property_map<arc_dissimilarity, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Arc, double> type;
    static const std::string name;
  };



  ////
  //// node_label
  ////
  struct node_label {};
  template <typename Graph>
  struct property_map<node_label, Graph> {
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, label_type> type;
    static const std::string name;
  };



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
    typedef boost::shared_ptr<AdjacencyList<T, U> > type;
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
  template <int N>
  class AdjacencyListBuilder {
  public:
    typedef boost::shared_ptr<NeighborhoodVisitorBase<N, label_type> > NeighborhoodVisitorPtr;
    typedef typename vigra::CoupledIteratorType<N, label_type, label_type>::type Iterator;
  private:
    vigra::MultiArrayView<N, unsigned> label_image_;
    NeighborhoodVisitorPtr neighborhood_accessor_;
    boost::shared_ptr<ListInserterBase<label_type, label_type> > list_inserter_;
  public:
    AdjacencyListBuilder();
    AdjacencyListBuilder(vigra::MultiArrayView<N, label_type> label_image,
                         NeighborhoodVisitorPtr neighborhood_accessor,
                         boost::shared_ptr<ListInserterBase<label_type, label_type> > list_inserter);
    void create_adjacency_list();
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
                       const typename AdjacencyListBuilder<N>::Iterator pixel) = 0;
  };


  ////
  //// NeighborhoodVisitor2D4Neighborhood
  ////
  template <typename T>
  class NeighborhoodVisitor2DRightLower : public NeighborhoodVisitorBase<2, T> {
  public:
    typedef boost::shared_ptr<PixelActorBase<T, T, typename AdjacencyList<T, T>::type > > PixelActorPtr;
  private:
    PixelActorPtr pixel_actor_;
    NeighborhoodVisitor2DRightLower();
  public:
    NeighborhoodVisitor2DRightLower(PixelActorPtr pixel_actor);
    virtual void visit(const vigra::MultiArrayView<2, T> image,
                       const AdjacencyListBuilder<2>::Iterator pixel);
  };


  ////
  //// PixelActorBase
  ////
  template <typename T, typename U, typename V>
  class PixelActorBase {
  private:
  public:
    virtual ~PixelActorBase() {}
    virtual void act(const T& pixel_value,
                     const U& comparison_value,
                     V& result) = 0;
  };


  ////
  //// PixelActorFindNeighbors
  ////
  template <typename T>
  class PixelActorFindNeighbors : public PixelActorBase<T, T, ListInserterBase<T, T> > {
  private:
  public:
    virtual void act(const T& pixel_value,
                     const T& comparison_value,
                     ListInserterBase<T, T>& result);
  };
  

  /* IMPLEMENTATIONS */
  

  ////
  //// node_neighbors
  ////  
  template <typename Graph>
  const std::string property_map<node_neighbors, Graph>::name = "node_neighbors";


  ////
  //// arc_dissimilarity
  ////
  template <typename Graph>
  const std::string property_map<arc_dissimilarity, Graph>::name = "arc_dissimilarity";


  ////
  //// node_label
  ////
  template <typename Graph>
  const std::string property_map<node_label, Graph>::name = "node_label";


  ////
  //// AdjacencyListBuilder
  ////
  template <int N>
  AdjacencyListBuilder<N>::AdjacencyListBuilder() :
    label_image_(),
    neighborhood_accessor_(),
    list_inserter_() {

  }


  template <int N>
  AdjacencyListBuilder<N>::AdjacencyListBuilder(vigra::MultiArrayView<N, unsigned> label_image,
                                                NeighborhoodVisitorPtr neighborhood_accessor,
                                                boost::shared_ptr<ListInserterBase<label_type, label_type> > list_inserter) :
    label_image_(label_image),
    neighborhood_accessor_(neighborhood_accessor),
    list_inserter_(list_inserter) {

  }


  template<int N>
  void AdjacencyListBuilder<N>::create_adjacency_list() {
    if (!label_image_.hasData()) {
      return;
    }
    // NeighborhoodVisitorBase<N, label_type>::PixelActorPtr pixel_actor(new PixelActorFindNeighbors<label_type>);
    Iterator start = createCoupledIterator(label_image_, label_image_);
    Iterator end = start.getEndIterator();
    for (Iterator it = start; it != end; ++it) {
      // do something
    }
  }


  ////
  //// ListInserterMap
  ////
  template <typename T, typename U>
  ListInserterMap<T, U>::ListInserterMap() :
    adjacency_list_(new typename AdjacencyList<T, T>::type) {

  }

  
  template <typename T, typename U>
  ListInserterMap<T, U>::ListInserterMap(typename AdjacencyListPtr<T, U>::type adjacency_list) :
    adjacency_list_(adjacency_list) {

  }


  template <typename T, typename U>
  void ListInserterMap<T, U>::add_to_list(const T& key,
                                          const U& value) {
    // (*adjacency_list_)[
  }


  ////
  //// NeighborhoodVisitor2DRightLower
  ////
  template <typename T>
  NeighborhoodVisitor2DRightLower<T>::NeighborhoodVisitor2DRightLower(PixelActorPtr pixel_actor) : 
    pixel_actor_(pixel_actor) {
    
  }



  template <typename T>
  void NeighborhoodVisitor2DRightLower<T>::visit(const vigra::MultiArrayView<2, T> image,
                                                 const AdjacencyListBuilder<2>::Iterator pixel) {
    return;
  }
  

  ////
  //// PixelActorFindNeighbors
  ////
  template <typename T>
  void PixelActorFindNeighbors<T>::act(const T& pixel_value,
                                       const T& comparison_value,
                                       ListInserterBase<T, T>& result) {
    return;
  }






}

#endif /* MULTI_HYPOTHESES_SEGMENTATION_H */
