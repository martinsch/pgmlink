#ifndef MULTI_HYPOTHESES_SEGMENTATION_H
#define MULTI_HYPOTHESES_SEGMENTATION_H

// stl
#include <vector>
#include <set>
#include <map>

// vigra
#include <vigra/multi_array.hxx>

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

  template <typename T, typename U>
  class PixelActorBase;

  typedef unsigned label_type;

  typedef boost::shared_ptr<MultiSegmenter> MultiSegmenterPtr;

  typedef boost::shared_ptr<std::vector<vigra::MultiArray<2, label_type> > > AssignmentListPtr;

  typedef std::map<unsigned, std::set<label_type> > AdjacencyList;

  typedef boost::shared_ptr<AdjacencyList > AdjacenyListPtr;

  



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
    typedef boost::shared_ptr<NeighborhoodVisitorBase<N, unsigned> > NeighborhoodVisitorPtr;
  private:
    vigra::MultiArrayView<N, unsigned> label_image_;
    NeighborhoodVisitorPtr neighborhood_accessor_;
  public:
    AdjacencyListBuilder();
    AdjacencyListBuilder(vigra::MultiArrayView<N, unsigned> label_image,
                         NeighborhoodVisitorPtr neighborhood_accessor);
    AdjacenyListPtr create_adjacency_list();
  };


  ////
  //// NeighborhoodVisitorBase
  ////
  template <int N, typename T>
  class NeighborhoodVisitorBase {
  public:
    typedef boost::shared_ptr<PixelActorBase<T, T> > PixelActorPtr;
  private:
    PixelActorPtr pixel_actor_;
  public:
    virtual void visit(const typename vigra::MultiArrayView<N, T>::difference_type& index) = 0;
  };


  ////
  //// NeighborhoodVisitor2D4Neighborhood
  ////
  class NeighborhoodVisitor2D4Neighborhood : public NeighborhoodVisitorBase<2, label_type> {
  private:
    NeighborhoodVisitor2D4Neighborhood();
  public:
    NeighborhoodVisitor2D4Neighborhood(PixelActorPtr pixel_actor);
    virtual void visit(const typename vigra::MultiArrayView<2, label_type>::difference_type& index);
  };


  ////
  //// PixelActorBase
  ////
  template <typename T, typename U>
  class PixelActorBase {
  private:
  public:
    virtual void act(const T& pixel_value, const U& comparison_value) = 0;
  };
  

  /* IMPLEMENTATIONS */


  ////
  //// NeighborhoodVisitorBase
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
  AdjacenyListPtr AdjacencyListBuilder<N>::create_adjacency_list() {
    AdjacenyListPtr adjacency_list(new AdjacencyList);
    if (!label_image_.hasData()) {
      return adjacency_list;
    }
    
    return adjacency_list;
  }
  
}

#endif /* MULTI_HYPOTHESES_SEGMENTATION_H */
