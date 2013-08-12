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
#include <vigra/algorithm.hxx> // argMax

// boost
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

// pgmlink
#include <pgmlink/log.h>
#include <pgmlink/region_graph.h>
#include <pgmlink/multi_hypotheses_segmentation.h>
#include <pgmlink/image_input.h>
#include <pgmlink/clustering.h>


namespace va = vigra::acc;

namespace pgmlink {

  template <typename T, int N>
  class MultiSegmentBuilder;

  template <typename T, int N>
  class MultiHypothesesGraphVectorBuilder;
  

  ////
  //// MultiSegmentBuilder
  ////
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
                        typename ImageWriterPtr<T, N>::type image_writer,
                        ClusteringBuilderPtr clustering_builder,
                        unsigned clusters_per_connected_component
                        );
  private:
    ImageSelectorPtr image_selector_;
    typename ImageRetrieverPtr<T, N>::type image_retriever_;
    typename ImageWriterPtr<T, N>::type image_writer_;
    ClusteringBuilderPtr clustering_builder_;
    std::vector<unsigned> clusters_per_connected_component_;
  };


  ////
  //// MultiHypothesesGraphVectorBuilder
  ////
  /* template <typename T, int N>
  class MultiHypothesesGraphVectorBuilder {
  public:
    boost::shared_ptr<std::vector<boost::shared_ptr<RegionGraph> > > build(const std::string& raw_directory,
                                                                           const std::string& label_directory);
    MultiHypothesesGraphVectorBuilder(
                                      ImageSelectorPtr image_selector,
                                      typename ImageRetrieverPtr<T, N>::type image_writer,
                                      );
  private:
    ClusteringBuilderPtr clustering_builder_;
    ImageSelectorPtr image_selector_;
    typename ImageRetrieverPtr<T, N>::type image_retriever_;
  }; */
                                
    
  
  /* IMPLEMENTATIONS */


  
  ////
  //// MultiSegmentBuilder
  ////
  template <typename T, int N>
  void MultiSegmentBuilder<T, N>::create_hypotheses(
                                              const std::string& input_directory,
                                              const std::string& output_directory
                                              ) {
    LOG(logDEBUG) << "MultiSegmenterBuilder<T, " << N << ">::create_hypotheses()";
    unsigned starting_index = 1;
    FilenameListPtr filenames = image_selector_->select(input_directory);
    for (FilenameList::iterator filename = filenames->begin();
         filename != filenames->end();
         ++filename) {
      image_retriever_->options_.set("filename", *filename);
      vigra::MultiArray<N, T> input_image = image_retriever_->retrieve();
      vigra::MultiArray<N, label_type> label_image(input_image.shape());
      LOG(logDEBUG1) << "MultiSegmenterBuilder<T, " << N << ">::create_hypotheses()"
                    << " -- input image shape: " << input_image.shape();
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
        const std::vector<CoordType >& coordinates =
          va::get<va::Coord<va::ValueList> >(accumulator, region);
        
        for (typename std::vector<CoordType >::const_iterator v_it = coordinates.begin();
             v_it != coordinates.end();
             ++v_it) {
          // for (typename CoordType::const_iterator c_it = coordinates.begin();
          // c_it != coordinates.end();
          // ++c_it) {
          // curr.push_back(*c_it);
          // }
          std::copy(v_it->begin(),
                    v_it->end(),
                    std::back_inserter(curr));
          if (N == 2) {
            curr.push_back(0);
          }
        }
      }
                                        
      std::string output_path = output_directory + "/" +
        boost::filesystem::path(*filename).filename().string();
      
      ConnectedComponentsToMultiSegments
        connected_components_to_multi_segments(
                                               components_coordinates,
                                               clusters_per_connected_component_,
                                               clustering_builder_,
                                               starting_index
                                               );
      vigra::MultiArray<4, label_type> output_image(vigra::Shape4(input_image.shape()[0],
                                                                  input_image.shape()[1],
                                                                  1,
                                                                  1)
                                                    );
      connected_components_to_multi_segments.to_images(output_image);
      starting_index = *vigra::argMax(output_image.begin(), output_image.end());
                                               
      image_writer_->options_.set("filename", output_path);
      image_writer_->write(output_image.bindAt(3, 0).bindAt(2, 0));
    }

  }


  template <typename T, int N>
  MultiSegmentBuilder<T, N>::MultiSegmentBuilder(
                                                 ImageSelectorPtr image_selector,
                                                 typename ImageRetrieverPtr<T, N>::type image_retriever,
                                                 typename ImageWriterPtr<T, N>::type image_writer,
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



  ////
  //// MultiHypothesesGraphVectorBuilder
  ////
  /* template <typename T, int N>
  boost::shared_ptr<std::vector<boost::shared_ptr<RegionGraph> > >
  MultiHypothesesGraphVectorBuilder<T, N>::build(const std::string& raw_directory,
                                                 const std::string& label_directory) {
    boost::shared_ptr<std::vector<boost::shared_ptr<RegionGraph> > > graphs;
    FilenameListPtr filenames = image_selector_->select(label_directory);
    std::sort(filenames.begin(), filenames.end());
    for (FilenameList::iterator filename = filenames->beign();
         filename != filenames->end();
         ++filename) {
      image_retriever_->options_.set("filename", *filename);
      vigra::MultiArray<N, T> input_image = image_retriever_->retrieve();
      AdjacencyGraphPtr graph(new AdjacencyGraph);
      ListInserterPtr<T, T>::type inserter(new ListInserterGraph<T, T>(graph));
      PixelActorPtr<T, T>::type actor(new PixelActorFindNeighbors<T>(inserter));
      NeighborhoodVisitorPtr<N, T>::type visitor(new NeighborhodVisitorAdjacentPlusOne<N, T>(actor));
      ParentConnectedComponentPtr<T, T>::type connected_component(new ParentConnectedComponentGraph<T, T>(graph));
      AdjacencyListBuilder<N, T> builder();
      builder.create_adjacency_list();
      graphs->push_back(graph);
    }
    return graphs;
  } */
}

#endif /* MULTI_HYPOTHESES_WORKFLOW_H */
