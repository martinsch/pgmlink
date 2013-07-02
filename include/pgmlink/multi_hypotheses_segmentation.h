#ifndef MULTI_HYPOTHESES_SEGMENTATION_H
#define MULTI_HYPOTHESES_SEGMENTATION_H

// stl
#include <vector>
#include <set>
#include <map>
#include <iostream>

// vigra
#include <vigra/multi_array.hxx>
#include <vigra/multi_iterator.hxx>
#include <vigra/multi_iterator_coupled.hxx>

// boost
#include <boost/shared_ptr.hpp>

// pgmlink
#include <pgmlink/clustering.h>
#include <pgmlink/traxels.h>


namespace pgmlink {
  class MultiSegmenter;
  
  class MultiSegmenterBuilder;

  class MultiSegmentContainer;

  class ConnectedComponentsToMultiSegments;

  template <int N>
  class AdjacencyListBuilder;
  
  template <int N, typename T>
  class NeighborhoodVisitorBase;

  template <typename T, typename U, typename V>
  class PixelActorBase;

  typedef unsigned label_type;

  typedef boost::shared_ptr<MultiSegmenter> MultiSegmenterPtr;

  typedef boost::shared_ptr<std::vector<vigra::MultiArray<2, label_type> > > AssignmentListPtr;

  typedef std::map<unsigned, std::set<label_type> > AdjacencyList;

  typedef boost::shared_ptr<AdjacencyList > AdjacencyListPtr;

  



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
  public:
    AdjacencyListBuilder();
    AdjacencyListBuilder(vigra::MultiArrayView<N, label_type> label_image,
                         NeighborhoodVisitorPtr neighborhood_accessor);
    AdjacencyListPtr create_adjacency_list();
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
    typedef boost::shared_ptr<PixelActorBase<T, T, AdjacencyList> > PixelActorPtr;
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
  class PixelActorFindNeighbors : public PixelActorBase<T, T, AdjacencyList> {
  private:
  public:
    virtual void act(const T& pixel_value,
                     const T& comparison_value,
                     AdjacencyList& result);
  };
  

  /* IMPLEMENTATIONS */


  ////
  //// AdjacencyListBuilder
  ////
  template <int N>
  AdjacencyListBuilder<N>::AdjacencyListBuilder() :
    label_image_(),
    neighborhood_accessor_() {

  }


  template <int N>
  AdjacencyListBuilder<N>::AdjacencyListBuilder(vigra::MultiArrayView<N, unsigned> label_image,
                                                NeighborhoodVisitorPtr neighborhood_accessor) :
    label_image_(label_image),
    neighborhood_accessor_(neighborhood_accessor) {

  }


  template<int N>
  AdjacencyListPtr AdjacencyListBuilder<N>::create_adjacency_list() {
    AdjacencyListPtr adjacency_list(new AdjacencyList);
    if (!label_image_.hasData()) {
      return adjacency_list;
    }
    // NeighborhoodVisitorBase<N, label_type>::PixelActorPtr pixel_actor(new PixelActorFindNeighbors<label_type>);
    Iterator start = createCoupledIterator(label_image_, label_image_);
    Iterator end = start.getEndIterator();
    for (Iterator it = start; it != end; ++it) {
      // do something
    }
    return adjacency_list;
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
                                       AdjacencyList& result) {
    return;
  }






}

#endif /* MULTI_HYPOTHESES_SEGMENTATION_H */
