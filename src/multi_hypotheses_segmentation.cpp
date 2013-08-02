// stl
#include <vector>
#include <numeric>

// vigra
#include <vigra/multi_array.hxx>

// boost
#include <boost/shared_ptr.hpp>

// armadillo
#include <armadillo>

// pgmlink
#include <pgmlink/multi_hypotheses_segmentation.h>
#include <pgmlink/clustering.h>
#include <pgmlink/traxels.h>
#include "pgmlink/log.h"

namespace pgmlink {
  ////
  //// MultiSegmenter
  ////
  MultiSegmenter::MultiSegmenter(const std::vector<unsigned>& n_clusters,
                                 ClusteringPtr clusterer) :
    n_clusters_(n_clusters),
    clusterer_(clusterer) {
    LOG(logDEBUG) << "MultiSegmenter::MultiSegmenter()";
  }


  vigra::MultiArray<2, label_type> MultiSegmenter::operator()(uint offset) {
    LOG(logDEBUG) << "MultiSegmenter::operator()(uint offset)";
    const arma::mat& data_arma_ = clusterer_->get_data_arma();
    unsigned n_samples = data_arma_.n_cols;
    unsigned n_layers = n_clusters_.size();
    vigra::MultiArray<2, label_type> res(vigra::Shape2(n_samples, n_layers), 0u);
    std::vector<unsigned>::const_iterator it = n_clusters_.begin();
    unsigned layer_index = 0;
    unsigned sample_index;
    for (; it != n_clusters_.end(); ++it, ++layer_index) {
      clusterer_->set_k_clusters(*it);
      clusterer_->operator()();
      for (sample_index = 0; sample_index < n_samples; ++sample_index) {
        res(sample_index, layer_index) = assign(data_arma_.col(sample_index)) + offset;
      }
      offset += *it;
    }
    return res;
  }


  label_type MultiSegmenter::assign(const arma::vec& sample) {
    LOG(logDEBUG) << "MultiSegmenter::assign(const arma::vec& sample)";
    return clusterer_->get_cluster_assignment(sample);
  }


  ////
  //// MultiSegmenterBuilder
  ////
  MultiSegmenterBuilder::MultiSegmenterBuilder(const std::vector<unsigned>& n_clusters,
                                               ClusteringBuilderPtr clustering_builder) :
    n_clusters_(n_clusters),
    clustering_builder_(clustering_builder) {
    LOG(logDEBUG) << "MultiSegmenterBuilder::MultiSegmenterBuilder(";
    LOG(logDEBUG) << "const std::vector<unsigned>& n_clusters,";
    LOG(logDEBUG) << "ClusteringBuilderPtr clustering_builder)";
  }


  MultiSegmenterPtr MultiSegmenterBuilder::build(const feature_array& data) {
    LOG(logDEBUG) << "MultiSegmenterBuidler::build(const feature_array& data)";
    return MultiSegmenterPtr(new MultiSegmenter(n_clusters_,
                                                clustering_builder_->build(data)
                                                )
                             );
  }


  ////
  //// MultiSegmentContainer
  ////
  MultiSegmentContainer::MultiSegmentContainer(vigra::MultiArrayView<2, label_type> assignments,
                                               const feature_array& coordinates) :
    assignments_(assignments), coordinates_(coordinates) {
    LOG(logDEBUG) << "MultiSegmentContainer::MultiSegmentContainer(";
    LOG(logDEBUG) << "vigra::MultiArrayView<2, label_type> assignments,";
    LOG(logDEBUG) << "const feature_array& coordinates)";
  }


  int MultiSegmentContainer::to_images(vigra::MultiArrayView<4, label_type> dest) {
    LOG(logDEBUG) << "MultiSegmentContainer::to_images(vigra::MultiArrayView<4, label_type> dest)";
    if (coordinates_.size()/3 != static_cast<unsigned>(assignments_.shape()[0])) {
      // number of samples in coordinates differ from number of samples
      // in assignments
      return 1;
    }
    
    if (assignments_.shape()[1] != dest.shape()[3]) {
      // number of segmentation layers does not fit destination image
      return 2;
    }
    
    unsigned n_layers = dest.shape()[3];
    unsigned n_samples = assignments_.shape()[0];
    unsigned sample_coord_index;
    unsigned coord1, coord2, coord3;
    for (unsigned layer = 0; layer < n_layers; ++layer) {
      for (unsigned sample = 0; sample < n_samples; ++sample) {
        sample_coord_index = 3*sample;
        coord1 = coordinates_[sample_coord_index++];
        coord2 = coordinates_[sample_coord_index++];
        coord3 = coordinates_[sample_coord_index];
        dest(coord1, coord2, coord3, layer) = assignments_(sample, layer);
      }
    }

    // regular return
    return 0;
  }


  ////
  //// ConnectedComponentsToMultiSegments
  ////
  ConnectedComponentsToMultiSegments::
  ConnectedComponentsToMultiSegments(const std::vector<feature_array >& components_coordinates,
                                     const std::vector<unsigned>& n_clusters,
                                     ClusteringBuilderPtr clustering_builder,
                                     unsigned starting_index) :
    components_coordinates_(components_coordinates),
    n_clusters_(n_clusters),
    clustering_builder_(clustering_builder),
    starting_index_(starting_index) {
    LOG(logDEBUG) << "ConnectedComponentsToMultiSegments::";
    LOG(logDEBUG) << "ConnectedComponentsToMultiSegments(";
    LOG(logDEBUG) << "const std::vector<feature_array >& components_coordinates,";
    LOG(logDEBUG) << "const std::vector<unsigned>& n_clusters,";
    LOG(logDEBUG) << "ClusteringBuilderPtr clustering_builder,";
    LOG(logDEBUG) << "unsigned starting_index)";
  }


  AssignmentListPtr ConnectedComponentsToMultiSegments::get_assignments() {
    LOG(logDEBUG) << "ConnectedComponentsToMultiSegments::get_assignments()";
    if (assignments_) {
      return assignments_;
    } else {
      assignments_ = AssignmentListPtr(new std::vector<vigra::MultiArray<2, unsigned> >);
      MultiSegmenterPtr segmenter;
      unsigned offset = starting_index_;
      unsigned offset_step = std::accumulate(n_clusters_.begin(), n_clusters_.end(), 0);
      std::vector<feature_array >::const_iterator cc_it = components_coordinates_.begin();
      for (; cc_it != components_coordinates_.end(); ++cc_it) {
        MultiSegmenterBuilder segmenter_builder(n_clusters_, clustering_builder_);
        segmenter = segmenter_builder.build(*cc_it);
        assignments_->push_back(segmenter->operator()(offset));
        offset += offset_step;
      }
      return assignments_;
    }
  }


  void ConnectedComponentsToMultiSegments::to_images(vigra::MultiArray<4, unsigned>& dest) {
    LOG(logDEBUG) << "ConnectedComponentsToMultiSegments::to_images(vigra::MultiArray<4, unsigned>& dest)";
    get_assignments();
    assert(components_coordinates_.size() == assignments_->size());
    std::vector<feature_array >::const_iterator cc_it = components_coordinates_.begin();
    std::vector<vigra::MultiArray<2, unsigned> >::iterator as_it = assignments_->begin();
    for (; cc_it != components_coordinates_.end(); ++cc_it, ++as_it) {
      MultiSegmentContainer segment_container(*as_it, *cc_it);
      segment_container.to_images(dest);
    }
  }


  ////
  //// RegionMergingGraph
  ////
  RegionMergingGraph::RegionMergingGraph() :
    adjacency_graph_(new AdjacencyGraph),
    maximum_merges_per_connected_component_(1u),
    maximum_merges_per_patch_(1u),
    starting_label_(0) {
    LOG(logDEBUG) << "RegionMergingGraph::RegionMergingGraph()";
    adjacency_graph_->add(node_merged_n_times());
  }


  RegionMergingGraph::RegionMergingGraph(AdjacencyGraphPtr adjacency_graph,
                                         unsigned maximum_merges_per_connected_component,
                                         unsigned maximum_merges_per_patch,
                                         label_type starting_label) :
    adjacency_graph_(adjacency_graph),
    maximum_merges_per_connected_component_(maximum_merges_per_connected_component),
    maximum_merges_per_patch_(maximum_merges_per_patch),
    starting_label_(starting_label) {
    LOG(logDEBUG) << "RegionMergingGraph::RegionMergingGraph(";
    LOG(logDEBUG) << "AdjacencyGraphPtr adjacency_graph,";
    LOG(logDEBUG) << "unsigned maximum_merges_per_connected_component,";
    LOG(logDEBUG) << "unsigned maximum_merges_per_patch,";
    LOG(logDEBUG) << "label_type starting_label)";
    adjacency_graph_->add(node_merged_n_times());
    // .add(connected_component_merged_n_times());
  }


  void RegionMergingGraph::merge() {
    //
    //
    // FIXME!!
    // SPLIT THIS HUGE BEAST UP INTO SMALLER FUNCTIONS
    // FIXME!!
    //
    //
    
    LOG(logDEBUG) << "void RegionMergingGraph::merge()";

    // get meta data containers for merging ready:
    // for each connected component the number of merges must not be larger than
    // the given maximum value
    std::map<label_type, unsigned> number_of_merges_per_connected_component;
    // store how often a region has been merged
    std::map<label_type, std::map<AdjacencyGraph::Region, int> > number_of_iterations_without_merge_global;

    AdjacencyGraph::ConnectedComponentMap& component_map = adjacency_graph_->get(node_connected_component());
    AdjacencyGraph::NeighborMap& neighbor_map = adjacency_graph_->get(node_neighbors());
    AdjacencyGraph::LabelMap& label_map = adjacency_graph_->get(node_label());

    // loop through components and merge regions until
    // number of merges reaches maximum for each component
    for (RegionGraph::ConnectedComponentMap::ValueIt component_it = component_map.beginValue();
           component_it != component_map.endValue();
           ++component_it) {
      if (*component_it == 0) {
        continue;
      }
      LOG(logDEBUG2) << "RegionMergingGraph::merge(): looping over components. Current component: " << *component_it;
      *component_it;
      unsigned& number_of_merges =
          number_of_merges_per_connected_component[*component_it];
      bool break_condition_component = true;
      std::map<AdjacencyGraph::Region, int>& number_of_iterations_without_merge =
          number_of_iterations_without_merge_global[*component_it];
      
      // initalize number_of_iterations_without_merge to zero
      LOG(logDEBUG3) << "RegionMergingGraph::merge() -- initializing number of iterations without merge to zero for each region";
      for (RegionGraph::ConnectedComponentMap::ItemIt item_it(component_map, *component_it);
           item_it != lemon::INVALID;
           ++item_it) {
        number_of_iterations_without_merge[item_it] = 0;
        LOG(logDEBUG4) << "RegionMergingGraph::merge() -- initialized number of iterations without merge to zero for region" << label_map[item_it];
      }

      // find region inside connected component that has been merged the least
      do {
        LOG(logDEBUG2) << "RegionMergingGraph::merge() -- merging regions until criteria are fulfilled for connected component "
                       << *component_it;
        std::map<AdjacencyGraph::Region, int>::iterator next_region_to_merge =
          number_of_iterations_without_merge.begin();
        for (std::map<AdjacencyGraph::Region, int>::iterator region_it =
               number_of_iterations_without_merge.begin();
             region_it != number_of_iterations_without_merge.end();
             ++region_it) {
          if (region_it->second < next_region_to_merge->second) {
            next_region_to_merge = region_it;
            
          }
          LOG(logDEBUG4) << "RegionMergingGraph::merge() -- current merging candidate: " << label_map[next_region_to_merge->first];
          // decrease (= increase absolute value) of count of
          // iterations without merge for that region
          region_it->second -= 1;
        }
        LOG(logDEBUG3) << "RegionMergingGraph::merge() -- next_region_to_merge: " << label_map[next_region_to_merge->first];
        std::set<AdjacencyGraph::Region>::iterator neighbor_to_merge =
          neighbor_map[next_region_to_merge->first].begin();
        int minimum = number_of_iterations_without_merge[*neighbor_to_merge];

        // find neighbor that has been merged the least of least merged region 
        for (std::set<AdjacencyGraph::Region>::iterator neighbor_it =
               neighbor_map[next_region_to_merge->first].begin();
             neighbor_it != neighbor_map[next_region_to_merge->first].end();
             ++neighbor_it) {
          if (minimum > number_of_iterations_without_merge[*neighbor_it]) {
            minimum = number_of_iterations_without_merge[*neighbor_it];
            neighbor_to_merge = neighbor_it;
          }
          number_of_iterations_without_merge[*neighbor_it] -= 1;
          LOG(logDEBUG4) << "RegionMergingGraph::merge() -- current neighbor merge candidate: " << label_map[*neighbor_to_merge];
        }
        LOG(logDEBUG3) << "RegionMergingGraph::merge() -- neighbor_to_merge: " << label_map[*neighbor_to_merge];

        // merge least merged region and least merged neighbor
        adjacency_graph_->merge_regions(next_region_to_merge->first, *neighbor_to_merge);

        // reset counter 
        next_region_to_merge->second = 0;
        number_of_iterations_without_merge[*neighbor_to_merge] = 0;
        number_of_merges += 1;

        // stop merging after connected component has reached maximum
        if (number_of_merges < maximum_merges_per_connected_component_) {
          break_condition_component = false;
        }
        if (break_condition_component) {
          LOG(logDEBUG2) << "RegionMergingGraph::merge() -- merged regions for connected component "
                         << *component_it;
          break;
        }
      } while(true);
    }
  }


}





