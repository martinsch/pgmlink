#ifndef MULTI_HYPOTHESES_WORKFLOW_H
#define MULTI_HYPOTHESES_WORKFLOW_H

// stl
#include <vector>
#include <string>
#include <algorithm> // std::copy
#include <iterator> // std::back_inserter

// vigra
#include <vigra/multi_array.hxx> // MultiArray(View)
#include <vigra/accumulator.hxx> // feature accumulators
#include <vigra/labelimage.hxx> // connected components

// boost
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

// pgmlink
#include <pgmlink/region_graph.h>
#include <pgmlink/multi_hypotheses_segmentation.h>
#include <pgmlink/image_input.h>
#include <pgmlink/clustering.h>


namespace va = vigra::acc;

namespace pgmlink {

  template <typename T, int N>
  class MultiSegmentBuilder;

  class MultiHypothesesGraphVectorBuilder;
  

  template <typename T, int N>
  class MultiSegmentBuilder {
  public:
    typedef typename vigra::CoupledIteratorType<2, T, label_type>::type Iterator;
    typedef typename Iterator::value_type Handle;
    typedef vigra::TinyVector<long, N> CoordType;
    void create_hypotheses(const std::string& input_directory,
                           const std::string& output_directory);
    MultiSegmentBuilder(
                        ImageSelectorPtr image_selector,
                        typename ImageRetrieverPtr<T, N>::type image_retriever,
                        ImageWriterPtr image_writer,
                        ClusteringBuilderPtr clustering_builder,
                        unsigned clusters_per_connected_component
                        );
  private:
    ImageSelectorPtr image_selector_;
    typename ImageRetrieverPtr<T, N>::type image_retriever_;
    ImageWriterPtr image_writer_;
    ClusteringBuilderPtr clustering_builder_;
    std::vector<unsigned> clusters_per_connected_component_;
  };


  /* class MultiHypothesesGraphVectorBuilder {
  public:
    boost::shared_ptr<RegionGraph> build(const std::string& directory);
    MultiHypothesesGraphVectorBuilder(
                                      ClusteringBuilderPtr clutersing_builder,
                                      ImageSelectorPtr image_selector,
                                      const PixelActorFactory& pixel_actor,
                                      const NeighborhoodVisitorFactory& neighborhood_visitor,
                                      const ParentConnectedComponentFactory& parent_connected_component
                                      );
  }; */
                                
    
  
  /* IMPLEMENTATIONS */


  template <typename T, int N>
  void MultiSegmentBuilder<T, N>::create_hypotheses(
                                              const std::string& input_directory,
                                              const std::string& output_directory
                                              ) {
    unsigned starting_index = 1;
    FilenameListPtr filenames = image_selector_->select(input_directory);
    for (FilenameList::iterator filename = filenames->begin();
         filename != filenames->end();
         ++filename) {
      image_retriever_->set_image(*filename);
      vigra::MultiArrayView<N, T> input_image = image_retriever_->retrieve();
      vigra::MultiArray<N, label_type> label_image(input_image.shape());
      int regions = vigra::labelImageWithBackground(srcImageRange(input_image),
                                                    destImage(label_image),
                                                    true,
                                                    0);
      va::AccumulatorChainArray<Handle,
                                va::Select<va::DataArg<1>, va::LabelArg<2>,
                                           va::Coord<va::ValueList>, va::Count
                                           >
                                > accumulator;

      Iterator start = vigra::createCoupledIterator(input_image, label_image);
      Iterator end = start.getEndIterator();
      accumulator.ignoreLabel(0);
      va::extractFeatures(start, end, accumulator);

      std::vector<feature_array > components_coordinates;
      for (int region = 1; region < regions; ++region) {
        components_coordinates.push_back(pgmlink::feature_array());
        feature_array& curr = *(components_coordinates.end()-1);
        const std::vector<CoordType>& coordinates =
          va::get<va::Coord<va::ValueList> >(accumulator, region);
        std::copy(coordinates.begin(),
                  coordinates.end(),
                  std::back_inserter(curr));
        if (N == 2) {
          curr.push_back(0);
        }
      }
                                        
      std::string output_path = output_directory +
        boost::filesystem::path(*filename).filename().string();
      
      ConnectedComponentsToMultiSegments
        connected_components_to_multi_segments(
                                               components_coordinates,
                                               clusters_per_connected_component_,
                                               clustering_builder_,
                                               starting_index
                                               );
                                               
      image_writer_->set_current_image(output_path);
    }

  }


  template <typename T, int N>
  MultiSegmentBuilder<T, N>::MultiSegmentBuilder(
                                                 ImageSelectorPtr image_selector,
                                                 typename ImageRetrieverPtr<T, N>::type image_retriever,
                                                 ImageWriterPtr image_writer,
                                                 ClusteringBuilderPtr clustering_builder,
                                                 unsigned clusters_per_connected_component
                                                 ) :
  image_selector_(image_selector),
    image_retriever_(image_retriever),
    image_writer_(image_writer),
    clustering_builder_(clustering_builder),
    clusters_per_connected_component_(1, clusters_per_connected_component) {
    // nothing to be done here
  }
                        
}

#endif /* MULTI_HYPOTHESES_WORKFLOW_H */
