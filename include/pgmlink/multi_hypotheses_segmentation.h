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

  class NeighborhoodVisitorBase;

  typedef unsigned label_type;

  typedef boost::shared_ptr<MultiSegmenter> MultiSegmenterPtr;

  typedef boost::shared_ptr<std::vector<vigra::MultiArray<2, label_type> > > AssignmentListPtr;

  typedef boost::shared_ptr<std::map<unsigned, std::set<label_type> > > AdjacenyListPtr;

  typedef boost::shared_ptr<NeighborhoodVisitorBase> NeighborhoodVisitorPtr;

  



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
  /* template <int N>
  class AdjacencyListBuilder {
  private:
    vigra::MultiArrayView<N, unsigned> label_image_;
    NeighborhoodVisitorPtr neighborhood_accessor;
  public:
    AdjacencyListBuilder();
    AdjacencyListBuilder(vigra::MultiArrayView<N, unsigned> label_image,
                         NeighborhoodVisitorPtr neighborhood_accessor);
    AdjacenyListPtr create_adjacency_list();
  }; */


  ////
  //// NeighborhoodVisitor
  ////
  

}

#endif /* MULTI_HYPOTHESES_SEGMENTATION_H */
