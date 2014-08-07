import optparse
import h5py
import sys
import numpy as np
sys.path.append('.')

import pgmlink

from oversegmentation import watershed3d as oversegment
from regionMerging import multilevelRegionMerging as merging
from traxelstoreGeneration import generateTraxelstore as store
from eventVector import eventVector

import vigra


def readHdf5dataset(fn, internal_path, t, ch, t_axis, ch_axis):
    with h5py.File(fn, 'r') as f:
       ds = f[internal_path]
       slicing = [slice(None),] * len(ds.shape)
       slicing[t_axis] = t
       slicing[ch_axis] = ch
       img = ds[tuple(slicing)]
    return img





if __name__ == "__main__":
   parser = optparse.OptionParser(usage="%prog [options] [predictionMap.h5]")  
   parser.add_option('--output', '-o', dest='out', default='./result.h5', type=str, help='output file (hdf5) [default=%default]')
   parser.add_option('--internalOut', dest='out_internal', default='exported_data', type=str, help='internal hdf5 path to output dataset [default=%default]')

   parser.add_option('--oversegmentationOnly', dest='seg_only', action="store_true", default=False, help='exit workflow after oversegmentation [default=%default]')

   parser.add_option('--outputOversegmentation', dest='out_os', default=None, type=str, help='output file (hdf5) for oversegmentation result [default=%default]')
   parser.add_option('--internalOversegmentation', dest='out_os_internal', default='exported_data', type=str, help='internal hdf5 path to oversegmentation dataset [default=%default]')
   
   parser.add_option('--outputConnectedComponents', dest='out_cc', default=None, type=str, help='output file (hdf5) for connected component result [default=%default]')
   parser.add_option('--internalConnectedComponents', dest='out_cc_internal', default='exported_data', type=str, help='internal hdf5 path to connected component dataset [default=%default]')
   
   # oversegmentation parameters:
   parser.add_option('--raw', '-r', dest='raw', default=None, type=str, help='h5 file of raw data; keep None if watershed should be computed on prediction map rather than on the raw data [default=%default]')
   parser.add_option('--internalPathRaw', dest='raw_internal', default='exported_data', type=str, help='internal hdf5 path to raw data [default=%default]')
   parser.add_option('--internalPathPredMap', dest='pred_internal', default='exported_data', type=str, help='internal hdf5 path to prediction map data [default=%default]')
   parser.add_option('--prediction-channel', dest='pred_ch', default=1, type=int, help='channel of prediction map containing foreground predictions [default=%default]')
   parser.add_option('--normalize-predictions', dest='normalize', action='store_true', help='normalize predictions to 0-255 before thresholding [default=%default]')
   parser.add_option('--sigma', '-s', dest='sigma', default=1.0, type=float, help='GaussianGradientMagnitude filter sigma on which the watershed is computed [default=%default]')
   parser.add_option('--background-threshold', dest='bg_threshold', default=128, type=float, help='cut off all values below this value as background [default=%default]')
   parser.add_option('--t-axis', dest='t_axis', default=0, type=int, help='index of time-axis [default=%default]')
   parser.add_option('--ch-axis', dest='ch_axis', default=-1, type=int, help='index of channel-axis [default=%default]')
   parser.add_option('--seedSigma', dest='seed_sigma', default=3.0, type=float, help='the sigma with which the gradient image should be smoothed to retrieve the seeds from, only applicable in the gradientSeededGradientMagnitude method [default=%default]')
   parser.add_option('--seedRadius', dest='seed_radius', default=2, type=int, help='disc radius of the seeds [default=%default]')
   parser.add_option('--sizeFilter', dest='size_filter', default=2, type=int, help='filter all objects smaller than this size [default=%default]')
   parser.add_option('--overseg-method', dest='method', default='dtSeededGradientMagnitude', type=str, help='Oversegmentation methods may be dtSeededGradientMagnitude or gradientSeededGradientMagnitude or dtSeededDt [default=%default]')

   # region merging parameters:
   parser.add_option('--percentile', dest='perc', default=0.8, type=float, help='percentile of edge weights to be merged in each level [default=%default]')
   parser.add_option('--progressive', dest='prog', default=1., type=float, help='progressive adaption of percentiles in each level [default=%default]')
   parser.add_option('--num-levels', dest='nlevels', default=4, type=int, help='number of region merging levels [default=%default]')
   
   # traxelstore generation:
   parser.add_option('--dumpTraxelstore', dest='dump', default=None, type=str, help='dump traxelstore to given file [default=%default]')
   parser.add_option('--loadTraxelstore', dest='load', default=None, type=str, help='load traxelstore from given file [default=%default]')
   parser.add_option('--detectionClassifier', dest='det_fn', default=None, type=str, help='hdf5 file for detection/count classifier [default=%default]')
   parser.add_option('--internalDetectionClassifier', dest='det_internal', default='/ObjectClassification/ClassifierForests/Forest0000', type=str, help='internal hdf5 path to the count probabilities [default=%default]')
   parser.add_option('--internalDetectionFeatureNames', dest='det_feature_names', default='/ObjectClassification/SelectedFeatures/Standard Object Features', type=str, help='internal hdf5 path to list of selected features [default=%default]')
   parser.add_option('--classifierFile', dest='rf_fn', default='', type=str, help='hdf5 file for region, transition, and division classifier [default=%default]')

   # speed-ups:
   parser.add_option('--tMin', dest='t_min', default=0, type=int, help='first timestep for segmentation and tracking [default=%default]')
   parser.add_option('--tMax', dest='t_max', default=None, type=int, help='end timestep (exclusive) for segmentation and tracking [default=%default]')
   parser.add_option('--serializeGraph', dest='serialize_to', default=None, type=str, help='serialize the hypotheses graph to file [default=%default]')
   parser.add_option('--deserializeGraph', dest='deserialize_from', default=None, type=str, help='deserialize the hypotheses graph from file [default=%default]')

   # tracking parameters
   parser.add_option('--det', dest='det', default=10, type=int, help='detection parameter [default=%default]')
   parser.add_option('--mov', dest='mov', default=10, type=int, help='move parameter [default=%default]')
   parser.add_option('--app', dest='app', default=100, type=int, help='appearance parameter [default=%default]')
   parser.add_option('--dis', dest='dis', default=100, type=int, help='disappearance parameter [default=%default]')
   parser.add_option('--count', dest='count', default=10, type=int, help='count parameter [default=%default]')
   parser.add_option('--div', dest='div', default=10, type=int, help='division parameter [default=%default]')
   parser.add_option('--onePerConflict', dest='onePerConflict', action="store_true", default=False, help='exactly one region per conflict set must be active [default=%default]')
   parser.add_option('--opp', dest='opp', default=10000, type=int, help='opportunity parameter [default=%default]')
   parser.add_option('--withoutConflictFactors', dest='without_conflict_factors', default=False, action='store_true', help='Do not use conflict factors instead of simple detection unaries.')
   parser.add_option('--withoutHierarchicalCount', dest='without_hierarchical_count', default=False, action='store_true', help='Do not use hierarchical count factor.')
   parser.add_option('--withoutCountingIncoming', dest='without_counting_incoming', default=False, action='store_true', help='Do not use counting incoming factor.')
   parser.add_option('--limitOutgoingArcs', dest='limit_outgoing_arcs', default=10, type=int, help='Limit number of outgoing arcs. The arcs will be pruned to the nearest conflict set in each connected component, if 0. [default=%default]')
   parser.add_option('--transitionAlpha', dest='transition_alpha', default=0, type=int, help='Squared distance tracking, if transition alpha is set to a value other than 0, this is alpha in the exponential decay then. In this case, the transition classifier is ignored. [default=%default]')
   parser.add_option('--maxDiv', dest='max_div', default=0, type=int, help='maximal division level [default=%default]')
   parser.add_option('--borderMargin', dest='border_margin', default=0, type=int, help='margin at border of field of view in which app/disapp costs are reduced linearly [default=%default]')
   parser.add_option('--epGap', dest='ep_gap', default=0.05, type=float, help='optimality gap for CPLEX [default=%default]')
   parser.add_option('--neighbors', dest='neighbors', default=1, type=int, help='number of neighbors considered for forward/backward linking [default=%default]')
   parser.add_option('--with-profiling', action='store_true', help='Activate profiling with yep', default=False)

   options, args = parser.parse_args()
   if len(args) != 1 and options.load is None and options.deserialize_from is None:
      parser.print_help()
      sys.exit(1)
            
   t_min = 0
   t_max = -1

   
   if options.with_profiling:
      import yep   
  
   if options.load is None and options.deserialize_from is None:

      with h5py.File(options.raw, 'r') as f:
         shape = list(f[options.raw_internal].shape)
      
      t_min = options.t_min
      t_max = options.t_max
      if t_max == None:
          t_max = shape[options.t_axis]
      
      if options.det_fn is None and not bool(options.seg_only):
         raise Exception, 'you need to give a path to the detection/count classifier'
     
      print 'preparing output files...'
      out_f = h5py.File(options.out, 'a')   
      out_shape = shape[:]
      out_shape[options.ch_axis] = options.nlevels + 1
      chunks = out_shape[:]
      chunks[options.t_axis] = 1
      chunks[options.ch_axis] = 1
      for idx, c in enumerate(chunks):
         if c > 64:
            chunks[idx] = 64
      try:
         del out_f[options.out_internal]
      except:
         pass
      out_ds = out_f.create_dataset(name=options.out_internal, shape=out_shape, chunks=tuple(chunks), compression=1, dtype=np.uint32)
 
      out_os_ds = None   
      if options.out_os is not None:
         out_os_f = h5py.File(options.out_os, 'a')
         try:
            del out_os_f[options.out_os_internal]
         except:
            pass
         out_os_ds = out_os_f.create_dataset(name=options.out_os_internal, shape=shape, chunks=tuple(chunks), compression=1, dtype=np.uint32)
 
      out_cc_ds = None   
      if options.out_cc is not None:
         out_cc_f = h5py.File(options.out_cc, 'a')
         try:
            del out_cc_f[options.out_cc_internal]
         except:
            pass
         out_cc_ds = out_cc_f.create_dataset(name=options.out_cc_internal, shape=shape, chunks=tuple(chunks), compression=1, dtype=np.uint32)
 
      connected_components = []
      conflict_sets = []
      region_maps = []
    
      for t in range(t_min, t_max):
         print
         print
         print 'oversegmentation at timestep =', t
         print '================================='

         predMap_at = readHdf5dataset(args[0], options.pred_internal, t, options.pred_ch, options.t_axis, options.ch_axis).squeeze()
         
         if bool(options.normalize):
            print '  normalizing prediction map...'
            predMap_at = oversegment.normalize(predMap_at)
         
         print '  thresholding prediction map (threshold =', options.bg_threshold, ')...'
         predMap_at_bin = oversegment.binarize(predMap_at, options.bg_threshold)

         if options.raw is None:
            img_at = predMap_at
         else:
            print '  read in raw dataset...'
            img_at = readHdf5dataset(options.raw, options.raw_internal, t, 0, options.t_axis, options.ch_axis).squeeze()

         assert predMap_at.shape == img_at.shape, 'the shapes of prediction image and raw image must match'
         
         print '  set background pixels (from thresholded prediction map) to zero...'
         img_at = oversegment.multiplyWith(img_at, predMap_at_bin, img_at.dtype)

         print '  running oversegmentation algorithm (sigma =', options.sigma, ') ...'
         res_with_bg_cluster, segmentImage_at = oversegment.do_oversegment(img_at, 
                                                           options.sigma, 
                                                           method=options.method, 
                                                           withOpening=False, # withOpening==True if an opening on on img_at should remove small objects
                                                           convert=False, # convertToRGB
                                                           seedSigma=options.seed_sigma,
                                                           seedDiscRadius=options.seed_radius,
                                                           sizeFilterFrom=options.size_filter
                                                           )
         
         print 
         print 'region merging at timestep =', t
         print '==============================='
      
         print '  building region adjacency graph...'
         regionGraph_at = merging.regionAdjacencyGraphFromSegmentImage(segmentImage_at)
         
         labelImages_at = []
         print '  drawing label image from region graph for level 0 (= segment image)...'
         alreadyAdded_at = []
         labelImages_at.append(merging.labelImageFromRegionGraph(regionGraph_at, segmentImage_at, alreadyAdded_at))

         for level in range(1, options.nlevels):
            perc = options.perc * pow(options.prog, level-1)
            print '  merging regions for level', level, 'with percentile=', perc, '...'
            regionGraph_at = merging.mergeNodesThreshold(regionGraph_at, percentile=perc)
      
            print '  drawing labelImage from region graph for level', level, '...'
            labelImages_at.append(merging.labelImageFromRegionGraph(regionGraph_at, segmentImage_at, alreadyAdded_at))

         print '  merging everything to connected components...'
         regionGraph_at = merging.mergeNodesThreshold(regionGraph_at, percentile=0.)
         print '  drawing label image from region graph for level', options.nlevels, '(= remaining connected components)'
         labelImages_at.append(merging.labelImageFromRegionGraph(regionGraph_at, segmentImage_at, alreadyAdded_at))

         connected_components_at = merging.getConnectedComponents(regionGraph_at)
         conflict_sets_at = merging.getConflictSets(regionGraph_at)
         region_maps_at = merging.getRegionToConnectedCompMap(regionGraph_at)
         assert len(region_maps_at) == max(region_maps_at.keys())         

         connected_components.append(connected_components_at)
         conflict_sets.append(conflict_sets_at)
         region_maps.append(region_maps_at)

         print 
         print 'writing results to file...'
         slicing = [slice(None),] * len(shape)
         slicing[options.t_axis] = t
         if out_os_ds is not None:
            out_os_ds[tuple(slicing)] = labelImages_at[0]

         if out_cc_ds is not None:
            out_cc_ds[tuple(slicing)] = merging.getConnectedComponentImage(region_maps_at, labelImages_at[0])

         for i in range(len(labelImages_at)):
            slicing[options.ch_axis] = slice(i,i+1)
            out_ds[tuple(slicing)] = labelImages_at[i]

     
      if out_os_ds is not None:
         out_os_f.close()
      if out_cc_ds is not None:
         out_cc_f.close()
      out_f.close()

      if bool(options.seg_only):
         print 'Oversegmentation done. Exit.'
         import sys; sys.exit(0)

      print
      print
      print 'generate traxelstore'
      print '===================='
      ts = pgmlink.MultiHypothesesTraxelStore()
      rf_detection = vigra.learning.RandomForest(options.det_fn, options.det_internal)
      with h5py.File(options.det_fn, 'r') as detection_file:
          detection_features = sorted(detection_file[options.det_feature_names].keys())

      # print ts
      classifier_h5 = options.rf_fn
      with h5py.File(classifier_h5, 'r') as f:
          features_move = f['Transitions']['Features']
          features_divisions = f['Divisions']['Features']
          features_detection = f['Regions']['Features']
          feature_names = set()
          for fn in list(features_move.values()) + list(features_divisions.values()) + list(features_detection.values()):
              fnames = list(fn.value)
              for ff in fnames:
                  #ff = ff.replace('object and ','')
                  for fff in ff.split(','):
                      feature_names.add(fff)
          if 'None' in feature_names:
               feature_names.remove('None')
          try:
              margin = tuple(f['margin'].value)
          except:
              margin = 0
      
      for t in range(t_min, t_max):
         print
         print '  timestep =', t
         
         store.generateTraxelStore_at(ts, t, conflict_sets[t], region_maps[t], options.out, options.out_internal, options.t_axis, options.ch_axis, options.nlevels, feature_names, rf_detection, detection_features, options.raw, options.raw_internal, margin=margin)
      
      if options.dump != None:
         import cPickle as pickle
         with open( options.dump, 'wb' ) as f_pickle:
             pickle.dump(ts, f_pickle )
   else: # options.load != None or options.deserialize_from != None
      if options.load != None:
          import cPickle as pickle
          with open( options.load, 'rb') as f_pickle:
              ts = pickle.load( f_pickle )

   # print ts
   classifier_h5 = options.rf_fn
   classifier_move = 'Transitions/rf'
   classifier_divisions = 'Divisions/rf'
   classifier_detection = 'Regions/rf'

   

   #rf_move = vigra.learning.RandomForest(classifier_h5, classifier_move)
   with h5py.File(classifier_h5, 'r') as f:
       features_move = f['Transitions']['Features']
       features_divisions = f['Divisions']['Features']
       features_detection = f['Regions']['Features']
   
       tracking_options = pgmlink.TrackingOptions()
       tracking_options.withDivisions(True)
       tracking_options.withDetectionVars(True)
       tracking_options.withConstraints(True)
       tracking_options.withMaximalConflictCliques(True)
       tracking_options.withConstantClassifiers(False)
       tracking_options.withClassifiers(True)
       tracking_options.withConstantClassifierFallback(True) # fall back to constant classifier if loading rf from h5 fails
       tracking_options.withClassifierCountPrecomputed(True) # count classifier already calculated in traxelstore generation
       tracking_options.forwardBackward(True) # multi hypotheses graph generation
       tracking_options.withConflictFactors(not(bool(options.without_conflict_factors)))
       tracking_options.set('det', options.det)
       tracking_options.set('div', options.div)
       tracking_options.set('app', options.app)
       tracking_options.set('dis', options.dis)
       tracking_options.set('forbidden', 0)
       tracking_options.set('max_div', options.max_div)
       tracking_options.set('distance', 1000000)
       tracking_options.set('gap', options.ep_gap)
       tracking_options.set('timeout', 1e+75)
       tracking_options.set('const_prob', 0.1)
       tracking_options.set('mov', options.mov)
       tracking_options.set('count', options.count)
       tracking_options.set('neighbors', options.neighbors)
       tracking_options.set('opportunity', options.opp)
       tracking_options.set('border_margin', options.border_margin)
       tracking_options.withOneActiveRegionPerComponent(options.onePerConflict) # by default, this should be false! if true, a hard constraint is added that at least one region in each conflict set must be active!
       ## field of view
       with h5py.File(options.raw, 'r') as f:
           shape = list(f[options.raw_internal].shape)
       if options.t_max is None or options.t_max == -1:
           t_max = shape[options.t_axis]
       else:
           t_max = options.t_max
       assert options.t_axis == 0, 'not implemented yet'
       assert options.ch_axis == len(shape)-1 or options.ch_axis == -1, 'not implemented yet'
       # assumes t,x,y,z,c
       if len(shape) - 2 < 3: # 2d
           z_to = 0
       elif len(shape) -2 == 3: # 3d
           z_to = shape[3]
       fov = pgmlink.FieldOfView(float(options.t_min), 0., 0., 0., t_max, float(shape[1]), float(shape[2]), float(z_to))
       tracking_options.fieldOfView(fov)
       ### SPEED-UPS
       if options.limit_outgoing_arcs > 0:
           tracking_options.limitOutgoingArcs( options.limit_outgoing_arcs )
       if options.transition_alpha > 0:
           tracking_options.withTransitionParameter( options.transition_alpha )
       tracking_options.withHierarchicalCountFactor( not(bool(options.without_hierarchical_count)) ) # count factor is hierarchical via hard constraints to avoid higher-order factors
       tracking_options.withCountingIncomingFactor( not(bool(options.without_counting_incoming)) ) # the incoming factor is always pairwise (counter variable is introduced if needed)

       # tracking_options.set('move', rf_move)
       tracking_options.set('classifier_file', classifier_h5)
       tracking_options.set('classifier_move', classifier_move)
       tracking_options.set('classifier_division', classifier_divisions)
       tracking_options.set('classifier_detection', classifier_detection)
       # tracking_options.set('classifier_count', '')

       if options.t_max is not None and options.t_max > -1:
           tracking_options.withRestrictedTimestepRange(options.t_min,options.t_max)

       for key in sorted(features_move.keys()):
           traxel_feats = features_move[key]
           for feat in sorted(traxel_feats):
               tracking_options.add('move', str(key), str(feat))
               print key, feat
               
       for key in sorted(features_divisions.keys()):
           traxel_feats = features_divisions[key]
           for feat in sorted(traxel_feats):
               tracking_options.add('division', str(key), str(feat))

       for key in sorted(features_detection.keys()):
           traxel_feats = features_detection[key]
           for feat in sorted(traxel_feats):
               tracking_options.add('detection', str(key), str(feat))



   if options.with_profiling:
      yep.start('out.prof')
   tracker = pgmlink.MultiHypothesesTracker(tracking_options)
   if options.deserialize_from is None:
       if options.serialize_to is None:
           events_vector = tracker.track(ts)
       else:
           events_vector = tracker.track(ts, str(options.serialize_to))
   else:
       if options.load != None:
           print 'WARNING: Ignoring the pickle traxelstore, since a hypotheses graph will be deserialized'
       print 'deserializing from', options.deserialize_from
       events_vector = tracker.track(str(options.deserialize_from))
   events = eventVector.get_events(events_vector)
   is_empty = True
   for t in sorted(events.keys()):
       print "t=%s" % t
       for key in events[t]:
           is_empty = False
           print "%d %s" % (len(events[t][key]), key)
   if options.out is not None:
       eventVector.write_events(events_vector, options.out)
       if is_empty == False:
          t_min = options.t_min
          t_max = options.t_max
          if t_max is None or t_max == -1:
               t_max = t_min + len(events)
          eventVector.project_active_file(options.out,
                                          TIME_AXIS = options.t_axis,
                                          LEVEL_AXIS = options.ch_axis,
                                          timestep_begin = t_min,
                                          timestep_end = t_max) # fix relabelling first!
   if options.with_profiling:
      yep.stop()
   


   


