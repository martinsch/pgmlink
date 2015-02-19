import pgmlink
import vigra
import numpy as np
import h5py

# for debugging
import inspect

def readHdf5dataset(fn, internal_path, t, ch, t_axis, ch_axis):
    with h5py.File(fn, 'r') as f:
       ds = f[internal_path]
       # if t is out of range, create a "fake" image containing only zeros
       # this is neccessary for the extraction of intersection features 
       if t == ds.shape[t_axis]:
           shape = list(ds.shape)
           del shape[t_axis]
           img = np.zeros(shape, dtype=ds.dtype)
       else:
           slicing = [slice(None),] * len(ds.shape)
           slicing[t_axis] = t
           slicing[ch_axis] = ch
           img = ds[tuple(slicing)]
    return img

def getVigraFeatures(ds_raw, ds, t_axis, ch_axis, features=['Count','RegionCenter','Mean','Variance'], margin=None):
   if t_axis < ch_axis:
      ch_axis -= 1 # since the t_axis is not included in ds
   if ds.shape[ch_axis] == 1:
      available_feats = vigra.analysis.supportedRegionFeatures(ds_raw.astype(np.float32), ds.astype(np.uint32))
      # parse vigra feature names
      available_feats_formatted = []
      for idx, ff in enumerate(available_feats):
          available_feats_formatted.append(ff.replace(' >', '>'))
      for idx,ff in enumerate(features):
        ff = ff.replace(' >','>')
        if ff not in available_feats_formatted:
            raise Exception, 'this feature is not available for this dataset!'
        
      f = getVigraFeaturesFromIlastik(ds_raw.astype(np.float32), ds.astype(np.uint32), features=available_feats, margin=margin)
      import sys; sys.exit(0)

      res = {}
      for key in f.keys():
         res[key.replace(' >','>')] = f[key]
      return res
   
   assert 'Count' in features, "'Count' feature must be computed in any case"

   feats_ch = []
   slicing = [ slice(None), ] * len(ds.shape)
   max_size = 1 # the first iteration will have a background label
   allUnique = np.unique(ds).tolist()

   slicing_raw = slicing[:]
   slicing_raw[ch_axis] = 0

   # extract features for each channel individually
   for ch in range(ds.shape[ch_axis]):
      slicing[ch_axis] = ch
      def op5ify(img):      
          if len(img.shape) == 2:
             img = img[np.newaxis,...,np.newaxis,np.newaxis]
          elif len(img.shape) == 3:
             img = img[np.newaxis,...,np.newaxis]
          elif len(ds_5d.shape) == 5:
             pass
          else:
             raise NotImplementedError, 'this shape is not handled yet'
          return img
      ds_5d = op5ify(ds[slicing].squeeze())
      ds_raw_5d = op5ify(ds_raw[slicing_raw].squeeze())

      # ilastik needs consecutive labels starting from 1
      uniqueLabels = np.unique(ds_5d).tolist()
      uniqueLabels.sort()
      if 0 in uniqueLabels:
          uniqueLabels.remove(0)
      count = 1
      relabelDict = {}
      for l in uniqueLabels:
          relabelDict[l] = count
          count += 1
      ds_5d = relabel(ds_5d, relabelDict)

      # parse available vigra features and get them from ilastik
      available_feats = vigra.analysis.supportedRegionFeatures(ds_raw_5d.astype(np.float32).squeeze(), ds_5d.astype(np.uint32).squeeze())
      available_feats_formatted = []
      for idx, ff in enumerate(available_feats):
          available_feats_formatted.append(ff.replace(' >', '>'))
      for idx,ff in enumerate(features):
        fff = ff.replace(' >','>').replace('object and ','').replace(' in neighborhood','')
        if fff not in available_feats_formatted :
            raise Exception, 'this feature is not available for this dataset: ' + str(ff)
        if 'in neighborhood' in ff:
            available_feats.append(ff)
      f = getVigraFeaturesFromIlastik(ds_raw_5d.astype(np.float32), ds_5d.astype(np.uint32), features=available_feats, margin=margin)
      assert f['Count'][0] == 0., ' the feature array must contain background'

      res = {}
      for key in f.keys():
          res[key.replace(' >','>')] = f[key]
      feats_ch.append(res)
      
      max_size += len(feats_ch[-1].values()[0]) - 1 # minus background

   # construct empty joint feature dictionary with max_size feature lists
   features = {}
   for feat_name in feats_ch[0].keys():
      features[feat_name] = [ 0 for i in range(max_size) ]

   # copy features into joint feature dictionary
   count = 1
   for feats in feats_ch:
      for row_idx in range(len(feats['Count'])):
         if feats['Count'][row_idx] == 0.:
            continue
         for feat_name in feats.keys():
            if not hasattr(feats[feat_name], '__iter__'):
                continue
            assert features[feat_name][count] == 0, "overwriting non-zero entry in joint feature dictionary"
            features[feat_name][count] = feats[feat_name][row_idx].flatten()
         count += 1
  
   assert len(features['Count']) == len(allUnique), 'inconsistent feature vector length'

   # filling background features and converting to numpy array
   assert features['Count'][1] != 0.
   for feat_name in features.keys():
      if hasattr(features[feat_name][1], '__iter__'):
         value = [ 0. for i in range(len(features[feat_name][1])) ]
         features[feat_name][0] = value
      features[feat_name] = np.array(features[feat_name])
         
   return features
         

def newTraxel(timestep, label, level):
   traxel = pgmlink.Traxel()
   traxel.Id = int(label)
   traxel.Timestep = int(timestep)
   traxel.Level = int(level)
   return traxel   

def getLevel(ds_at, region_id):
   return -1

def replace_missing(a):
   idxs = np.where(np.isnan(a) + np.isinf(a))
   if len(idxs[0]) == 0:
      return
   rows = list(set(idxs[0].flat))
   if len(idxs) == 2:
      cols = list(set(idxs[1].flat))
      idx = (rows, cols)
   else: 
      idx = rows
   a[idx] = 0
   return 

def addFeatureArray(traxel, features, label, feature_names):
   for name in feature_names:
      feats = features[name.replace(' >','>')][label]
      if isinstance(feats, np.ndarray):
         replace_missing(feats)
      feats = feats.tolist()
      if not isinstance(feats,list):
         feats = [feats]
      traxel.add_feature_array(str(name.replace(' >','>')), len(feats))
      for i, v in enumerate(feats):
         traxel.set_feature_value(str(name.replace(' >','>')), i, float(v))

   # add com feature   
   feats = features['RegionCenter'][label]
   if len(feats) == 2:
       feats = np.append(feats, np.zeros(1))
   traxel.add_feature_array("com", len(feats))
   for i, v in enumerate(feats):
      traxel.set_feature_value("com", i, float(v))

   return traxel

def getDetProbs(rf_det, features, region_id, conflictsMap):
   features = np.require(features, dtype=np.float32)
   probs = list(rf_det.predictProbabilities(features).flat)
   lenResult = len(conflictsMap[region_id]) + 1
   if lenResult < len(probs):
      result = probs[:lenResult-1]
      result.append(sum(probs[lenResult-1:]))      
   elif lenResult > len(probs):
      result = probs[:]
      for i in range(len(probs), lenResult):
         result.append(probs[-1])
      '''# renormalize
      const = sum(result)
      assert const != 0
      for i, v in enumerate(result):
         result[i] = float(v) / const'''
   else:
      result = probs
   assert len(result) == lenResult
   if sum(result) == 0:
      for i in range(len(result)):
         result[i] = 1./len(result)
   return result

def generateTraxelStore_at(ts, timestep, conflictsMap, regionToConnectedCompMap, fn, internal_path, t_axis, ch_axis, num_levels, feat_names, rf_det, feature_names_det, fn_raw, internal_raw, withDetProbs=True,margin=(2,2,2)):
   print '    reading in h5 dataset...'
   features_formatted = set([])
   for feature_name in feat_names:
       feature_name = feature_name.replace(' >', '>')
       features_formatted.add(feature_name)
   for feature_name in feature_names_det:
       feature_name = feature_name.replace(' >', '>')
       features_formatted.add(feature_name)
   features_formatted = list(features_formatted)
   ds_at = readHdf5dataset(fn, internal_path, timestep, slice(None), t_axis, ch_axis)
   ds_at_plus_one = readHdf5dataset(fn, internal_path, timestep+1, slice(None), t_axis, ch_axis)
   ds_raw_at = readHdf5dataset(fn_raw, internal_raw, timestep, slice(None), t_axis, ch_axis)
   print '    computing vigra features...'
   assert len(np.unique(ds_at)) == np.unique(ds_at)[-1]+1
   feats_at = getVigraFeatures(ds_raw_at, ds_at, t_axis, ch_axis, features_formatted, margin=margin)
   assert len(feats_at['Count']) - 1 == np.max(ds_at)
   print '    extracting intersections'
   ch_axis_sub = ch_axis
   if t_axis >= 0 and ch_axis > 0 and t_axis < ch_axis:
       ch_axis_sub = ch_axis - 1
   n_channels = ds_at.shape[ch_axis_sub]
   intersects = pgmlink.IntersectsDictionary()
   slicing1 = [slice(None)]*len(ds_at.shape)
   slicing2 = [slice(None)]*len(ds_at.shape)
   for channel1 in xrange(n_channels):
       slicing1[ch_axis_sub] = channel1
       for channel2 in xrange(n_channels):
           slicing2[ch_axis_sub] = channel2
           intersects.getIntersects(ds_at[tuple(slicing1)], ds_at_plus_one[tuple(slicing2)])

   print '    adding traxels...'
   for region_id in regionToConnectedCompMap.keys():
      trax = newTraxel(timestep, region_id, getLevel(ds_at, region_id))
      trax.Component = regionToConnectedCompMap[region_id]
      addFeatureArray(trax, feats_at, trax.Id, features_formatted)
      intersects.addToTraxel(trax, trax.Id)

      # add detProb for connected components only
      if trax.Component == region_id and withDetProbs:
         features = np.array([], dtype = np.float32)
         for feature_name in feature_names_det:
             feature_name = feature_name.replace(' >','>')
             features = np.append(features, feats_at[feature_name][region_id])
         probs = getDetProbs(rf_det, features[np.newaxis,...], region_id, conflictsMap)
         trax.add_feature_array("count_prediction", len(probs))

      ts.add(trax, regionToConnectedCompMap[region_id])

   print '    adding conflict maps'
   ts.addConflictMap(timestep, conflictsMap)

   return ts


def getVigraFeaturesFromIlastik(ds_raw, ds, features, ndim=None, margin=(2,2,2)):
    r = ds_raw.squeeze()
    l = ds.squeeze()    
    if ndim!=None:
        if len(r.shape) != ndim:
            print 'WARNING: this image does not have the right shape, adding newaxis'
        for i in range(ndim-len(r.shape)):
            r = r[...,np.newaxis]
            l = l[...,np.newaxis]
    from lazyflow.graph import Graph
    from lazyflow.operators import OpLabelImage
    from ilastik.applets.objectExtraction.opObjectExtraction import OpAdaptTimeListRoi, OpRegionFeatures
    from ilastik.plugins import pluginManager

    r = r.view(vigra.VigraArray)
    l = l.view(vigra.VigraArray)
    def setupDataset(img):
        # squeeze the channel
        img = img.squeeze()
        if len(img.shape) == 2:
            img = img[np.newaxis,...,np.newaxis,np.newaxis]
            img.axistags = vigra.defaultAxistags('txyzc')
        elif len(img.shape) == 3:
            img = img[np.newaxis,...,np.newaxis]
            img.axistags = vigra.defaultAxistags('txyzc')
        else:
            raise NotImplementedError
        return img
    r = setupDataset(r)
    l = setupDataset(l)
    g = Graph()
    op = OpRegionFeatures(graph=g)
    op.LabelImage.setValue(l)
    op.RawImage.setValue(r)
    vigra_name = "Standard Object Features"
    feat_names = { vigra_name: {} }
    for fname in features:
        #fname = fname.replace('>>>','> > >').replace('>>', '> >')
        fname = fname.replace(' >','>')
        fname_split = fname.split(',')
        assert len(fname_split) <= 2
        for ff in fname_split:
            ff = ff.replace('object and ', '')    
            if "neighborhood" in ff and margin != None:
                feat_names[vigra_name][ff] = {"margin": margin}        
            else:
                feat_names[vigra_name][ff] = {}
    op.Features.setValue(feat_names)
    op.Output.fixed = False
    #opAdapt = OpAdaptTimeListRoi(graph=op.graph)
    #opAdapt.Input.connect(op.Output)
    #feats = opAdapt.Output([0]).wait()
    feats = op.Output([]).wait()
    result = feats[0][vigra_name]
    assert len(np.unique(l)) == len(result.values()[0])
    return result

def relabel( volume, replace, defaultValue = 0 ):
   mp = np.arange(0,np.amax(volume)+1, dtype=volume.dtype)
   mp[1:] = defaultValue
   mp[replace.keys()] = replace.values()
   return mp[volume.astype(np.uint32)]

if __name__ == "__main__":
   import os
   conflictsMap_at = pgmlink.ConflictSetsMap()
   ts = pgmlink.MultiHypothesesTraxelStore()
   
   fn = os.path.abspath('../testData/mitocheck_overseg-merged.h5')
   fn_detProb = os.path.abspath('../testData/mitocheck_count-classifier.ilp')
   fn_raw = os.path.abspath('../testData/mitocheck_raw.h5')
   internal_raw = 'exported_data'
   ch_axis = 3
   t_axis = 2

   feature_names = ['RegionCenter',
                       'Count',
                       'Variance',
                       'Sum',
                       'Mean',
                       'RegionRadii',
                       'Central<PowerSum<2> >',
                       'Central<PowerSum<3> >',
                       'Central<PowerSum<4> >',
                       'Histogram',
                       'Kurtosis',
                       'Maximum',
                       'Minimum',
                       'Quantiles',
                       'RegionAxes',
                       'Skewness',
                       'Weighted<PowerSum<0> >']

   with h5py.File(fn,'r') as f:
      num_levels = f['exported_data'].shape[ch_axis]

   regionToConnectedCompMap = {}
   for i in range(595):
      regionToConnectedCompMap[i] = i

   rf_det = vigra.learning.RandomForest(fn_detProb, '/ObjectClassification/ClassifierForests/Forest0000')
   with h5py.File(fn_detProb, 'r') as detection_file:
       feature_names_det = sorted(detection_file['/ObjectClassification/SelectedFeatures/Standard Object Features'].keys())
   generateTraxelStore_at(ts,
                          0,
                          conflictsMap_at,
                          regionToConnectedCompMap,
                          fn,
                          'exported_data',
                          t_axis,
                          ch_axis,
                          num_levels,
                          feature_names,
                          rf_det,
                          feature_names_det,
                          fn_raw,
                          internal_raw,
                          False)
   trax = ts.get(0,268,268)
   indices = trax.features['intersect_indices']
   intersects = trax.features['intersects']
   assert 266 in indices
   assert 277 in indices
   assert 582 in indices
   for idx in range(indices.size()):
       print "Intersect of %d and %d: %d" % (268, indices[idx], intersects[idx])
