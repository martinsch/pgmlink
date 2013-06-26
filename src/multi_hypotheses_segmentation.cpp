// stl
#include <vector>

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
    clusterer_(clusterer) {}


  vigra::MultiArray<2, unsigned> MultiSegmenter::operator()(uint offset) {
    const arma::mat& data_arma_ = clusterer_->get_data_arma();
    unsigned n_samples = data_arma_.n_cols;
    unsigned n_layers = n_clusters_.size();
    vigra::MultiArray<2, unsigned> res(vigra::Shape2(n_samples, n_layers), 0u);
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


  unsigned MultiSegmenter::assign(const arma::vec& sample) {
    return clusterer_->get_cluster_assignment(sample);
  }


  ////
  //// MultiSegmenterBuilder
  ////
  MultiSegmenterBuilder::MultiSegmenterBuilder(const std::vector<unsigned>& n_clusters,
                                               ClusteringBuilderPtr clustering_builder) :
    n_clusters_(n_clusters),
    clustering_builder_(clustering_builder) {}


  MultiSegmenterPtr MultiSegmenterBuilder::build(const feature_array& data) {
    return MultiSegmenterPtr(new MultiSegmenter(n_clusters_,
                                                clustering_builder_->build(data)
                                                )
                             );
  }
}
