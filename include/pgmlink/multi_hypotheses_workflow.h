#ifndef MULTI_HYPOTHESES_WORKFLOW_H
#define MULTI_HYPOTHESES_WORKFLOW_H

// stl
#include <vector>
#include <string>
#include <algorithm> // std::copy, transform
#include <functional> //std::plus
#include <iterator> // std::back_inserter

// vigra
#include <vigra/multi_array.hxx> // MultiArray(View)
#include <vigra/accumulator.hxx> // feature accumulators
#include <vigra/labelimage.hxx> // connected components
#include <vigra/algorithm.hxx> // argMax
#include <vigra/accumulator.hxx> // feature accumulators

// boost
#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

// pgmlink
#include <pgmlink/log.h>
#include <pgmlink/region_graph.h>
#include <pgmlink/multi_hypotheses_segmentation.h>
#include <pgmlink/image_input.h>
#include <pgmlink/clustering.h>
#include <pgmlink/hypotheses.h>
#include <pgmlink/traxels.h>
#include <pgmlink/multi_hypotheses_graph.h>



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
                           const std::string& output_directory,
                           const std::string& label_directory);
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
  template <typename T, int N>
  class MultiHypothesesGraphVectorBuilder {
  public:
    typedef typename vigra::CoupledIteratorType<N, T, T>::type CoupledIterator;
    typedef typename CoupledIterator::value_type Handle;
    RegionGraphVectorPtr build(const std::string& raw_directory,
                                                                           const std::string& label_directory);
    MultiHypothesesGraphVectorBuilder(
                                      ImageSelectorPtr image_selector,
                                      typename ImageRetrieverPtr<T, N>::type image_retriever,
                                      unsigned maximum_merges_per_connected_component,
                                      unsigned maximum_merges_per_patch
                                      );
  private:
    ImageSelectorPtr image_selector_;
    typename ImageRetrieverPtr<T, N>::type image_retriever_;
    unsigned maximum_merges_per_connected_component_;
    unsigned maximum_merges_per_patch_;
  };
                                
    
  
  /* IMPLEMENTATIONS */


  
  ////
  //// MultiSegmentBuilder
  ////
  template <typename T, int N>
  void MultiSegmentBuilder<T, N>::create_hypotheses(
                                              const std::string& input_directory,
                                              const std::string& output_directory,
                                              const std::string& label_directory
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

      // make sure cluster regions have index at least greater than the
      // largest original region index
      starting_index = regions;
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


      std::string label_image_path = label_directory + "/" +
        boost::filesystem::path(*filename).filename().string();
      image_writer_->options_.set("filename", label_image_path);
      image_writer_->write(label_image);
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
  template <typename T, int N>
  RegionGraphVectorPtr
  MultiHypothesesGraphVectorBuilder<T, N>::build(const std::string& raw_directory,
                                                 const std::string& label_directory) {
    LOG(logDEBUG) << "MultiHypothesesGraphVectorBuilder<T, " << N << ">::build()";
    RegionGraphVectorPtr
      graphs(new std::vector<boost::shared_ptr<RegionGraph> >);
    FilenameListPtr filenames = image_selector_->select(label_directory);
    FilenameListPtr filenames_raw = image_selector_->select(raw_directory);
    std::sort(filenames->begin(), filenames->end());
    std::sort(filenames_raw->begin(), filenames_raw->end());
    FilenameList::iterator filename_raw = filenames_raw->begin();
    unsigned timestep = 0;
    for (FilenameList::iterator filename = filenames->begin();
         filename != filenames->end();
         ++filename, ++filename_raw, ++timestep) {
      image_retriever_->options_.set("filename", *filename);
      vigra::MultiArray<N, T> input_image = image_retriever_->retrieve();
      image_retriever_->options_.set("filename", *filename_raw);
      vigra::MultiArray<N, T> label_image = image_retriever_->retrieve();
      AdjacencyGraphPtr graph(new AdjacencyGraph);
      typename ListInserterPtr<T, T>::type inserter(new ListInserterGraph<T, T>(graph));
      typename PixelActorPtr<T, T>::type actor(new PixelActorFindNeighbors<T>(inserter));
      typename NeighborhoodVisitorPtr<N, T>::type visitor(new NeighborhoodVisitorAdjacentPlusOne<N, T>(actor));
      typename ParentConnectedComponentPtr<T, T>::type connected_component(new ParentConnectedComponentGraph<T, T>(graph));
      AdjacencyListBuilder<N, T> builder(input_image,
                                         label_image,
                                         visitor,
                                         connected_component);
      builder.create_adjacency_list();
      va::AccumulatorChainArray<Handle,
                                       va::Select<va::DataArg<1>, va::LabelArg<2>,
                                                  va::Mean, va::Variance, va::Coord<va::Mean>,
                                                  va::Count
                                                  >
                                       >
        accumulator_connected_components, accumulator_split_regions, accumulator;
      CoupledIterator start = vigra::createCoupledIterator(label_image, label_image);
      CoupledIterator end = start.getEndIterator();
      accumulator_connected_components.ignoreLabel(0);
      va::extractFeatures(start, end, accumulator_connected_components);
      start = vigra::createCoupledIterator(input_image, input_image);
      end = start.getEndIterator();
      accumulator_split_regions.ignoreLabel(0);
      va::extractFeatures(start, end, accumulator_split_regions);

      LOG(logDEBUG3) << "MultiHypothesesGraphVectorBuilder<T, " << N << ">::build()"
                     << " about to split merge accumulators with "
                     << accumulator_split_regions.maxRegionLabel()
                     << " and " << accumulator_connected_components.maxRegionLabel() << "labels";

      accumulator.setMaxRegionLabel(accumulator_split_regions.maxRegionLabel());
      accumulator_connected_components.setMaxRegionLabel(accumulator_split_regions.maxRegionLabel());
      accumulator += accumulator_connected_components;
      accumulator += accumulator_split_regions;

      graph->add(node_traxel());
      RegionGraph::TraxelMap& traxel_map = graph->get(node_traxel());
      RegionGraph::LabelMap& label_map = graph->get(node_label());
      RegionGraph::LevelMap& level_map = graph->get(node_level());
      RegionGraph::ContainingMap& containing_map = graph->get(node_contains());
      for (int i = 1; i <= accumulator.maxRegionLabel(); ++i) {
        Traxel trax(i, timestep);
        trax.features["com"] =
          feature_array(va::get<va::Coord<va::Mean> >(accumulator, i).begin(),
                        va::get<va::Coord<va::Mean> >(accumulator, i).end()
                        );
        if (trax.features["com"].size() == 2) {
          trax.features["com"].push_back(0.);
        }
        trax.features["size"] =
          feature_array(1, va::get<va::Count>(accumulator, i));
        const RegionGraph::Node& node = label_map(i);
        traxel_map.set(node, trax);
        LOG(logDEBUG3) << "MultiHypothesesGraphVectorBuilder<T, N>::build() -- " << trax
                       << " " << trax.features["com"][0] << "," << trax.features["com"][1]
                       << "," << trax.features["com"][2] << "; " << trax.features["size"][0];
      }
      
      RegionGraph::LabelMap::ValueIt value_iterator = label_map.beginValue();
      for (; value_iterator != label_map.endValue(); ++value_iterator) {
        if (*value_iterator <= accumulator.maxRegionLabel()) {
          continue;
        }
        const RegionGraph::Node& node = label_map(*value_iterator);
        const std::set<RegionGraph::Node>& children = containing_map[node];
        feature_array com(N, 0);
        feature_array size(1, 0);
        for (std::set<RegionGraph::Node>::const_iterator children_it = children.begin();
             children_it != children.end();
             ++children_it) {
          if (level_map[*children_it] > 0) {
            Traxel child_trax = traxel_map[*children_it];
            size[0] += child_trax.features["size"][0];
            const feature_array& child_com = child_trax.features["com"];
            std::transform(child_com.begin(), child_com.end(),
                           com.begin(), com.begin(),
                           std::plus<feature_type>()
                           );
          }
        }
        
        std::transform(com.begin(), com.end(), com.begin(),
                       std::bind2nd(std::divides<feature_type>(),
                                    size[0]
                                    )
                       );
        Traxel trax(*value_iterator, timestep);
        LOG(logDEBUG3) << "MultiHypothesesGraphVectorBuilder<T, N>::build() -- " << trax
                     << " " << com[0] << "," << com[1] << "," << com[2] << "; " << size[0];
        trax.features["size"] = size;
        trax.features["com"] = com;
        LOG(logDEBUG3) << "MultiHypothesesGraphVectorBuilder<T, N>::build() -- merged region size: " << size[0];
        
        
        traxel_map.set(node, trax);
      }


      graphs->push_back(graph);
    }



    return graphs;
  }


  template <typename T, int N>
  MultiHypothesesGraphVectorBuilder<T, N>::MultiHypothesesGraphVectorBuilder(ImageSelectorPtr image_selector,
                                                                             typename ImageRetrieverPtr<T, N>::type image_retriever,
                                                                             unsigned maximum_merges_per_connected_component,
                                                                             unsigned maximum_merges_per_patch
                                                                             ) :
    image_selector_(image_selector),
    image_retriever_(image_retriever),
    maximum_merges_per_connected_component_(maximum_merges_per_connected_component),
    maximum_merges_per_patch_(maximum_merges_per_patch) {
    // nothing to be done here
  }
}

#endif /* MULTI_HYPOTHESES_WORKFLOW_H */

