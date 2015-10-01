import numpy as np
import h5py
import vigra
import pgmlink
import datetime
import optparse

def getAnnotations(fn, path='/ManualTransitions/Labels/0', ):
   segment2Label = {}
   label2Segment = {}

   with h5py.File(fn, 'r') as f:
      labels = f[path]
      for t in labels.keys():
         labels_at = labels[t]
         t = int(t)
         if t not in segment2Label.keys():
            segment2Label[t] = {}
         if t not in label2Segment.keys():
            label2Segment[t] = {}

         for segment in labels_at.keys():
            segment2Label[t][int(float(segment))] = []
            for label in np.array(labels_at[segment]):
               segment2Label[t][int(float(segment))].append(int(float(label)))
               if int(float(label)) not in label2Segment[t].keys():
                  label2Segment[t][int(float(label))] = []
               label2Segment[t][int(float(label))].append(int(float(segment)))

   return segment2Label, label2Segment


def getDivisionAnnotations(fn, path='/ManualDivisions/Divisions/0'):
   positiveDivisions = []
   negativeDivisions = []
   with h5py.File(fn, 'r') as f:
      divisions = np.array(f[path])
      for row in divisions:
         if row[0] < 0:
            negativeDivisions.append(row[:-1].tolist())
         elif row[0] > 0:
            positiveDivisions.append(row[:-1].tolist())
         else:
            raise Exception, "shouldn't get here"

   return positiveDivisions, negativeDivisions


def getDivisionTripletsIlastik(fn, path='DivisionDetection/LabelInputs/0', labelImagePath='/TrackingFeatureExtraction/LabelImage/0', 
            featsPath='/TrackingFeatureExtraction/RegionFeaturesVigra/0/[[%d], [%d]]/Standard Object Features/', size_filter=3,
            templateSize=50):
    positiveDivisions = {}
    negativeDivisions = {}
    with h5py.File(fn, 'r') as f:
        labelImg_name = f[labelImagePath].keys()[0]
        # replace first 0 and first 1 by %d
        labelImg_name = labelImg_name.replace('0','%d', 1).replace('1', '%d', 1)        
        for t in f[path].keys():
            t = int(t)
            positiveDivisions[t] = []
            negativeDivisions[t] = []
            labels_at = f[path][str(t)].value
            for idx, val in enumerate(labels_at):
                val = int(val.flatten())
                if val == 0: # not annotated
                    continue
                li_cur_name = labelImg_name % (t, t+1)
                labelImg_cur = np.array(f[labelImagePath][li_cur_name]).squeeze()
                if t+1 < len(f[path].keys()):
                    li_next_name = labelImg_name % (t+1, t+2)
                    labelImg_next = np.array(f[labelImagePath][li_next_name]).squeeze()
                else:
                    #labelImg_next = np.zeros(labelImg_cur.shape)
                    continue
                
                feats_cur = f[featsPath % (t,t+1)]
                feats_next = f[featsPath % (t+1,t+2)]
                coords = feats_cur['RegionCenter'][idx]
                roi = len(coords) * [slice(None),]
                for i in range(len(coords)):
                    roi[i] = slice(max(0, coords[i]-templateSize/2), min(coords[i]+templateSize/2, labelImg_next.shape[i]))
                nearestLabels_next = np.unique(labelImg_next[roi]).tolist()
                if 0 in nearestLabels_next:
                    nearestLabels_next.remove(0)

                if len(nearestLabels_next) < 2:
                    print 'WARNING: not enough labels in neighborhood, increase the template size!'
                    continue
                
                nearestLabels_next_filtered = []
                distances = []
                for l in nearestLabels_next:
                    if feats_next['Count'][l] > size_filter:
                        nearestLabels_next_filtered.append(l)
                        distances.append(np.linalg.norm(feats_next['RegionCenter'][l] - coords))

                if len(nearestLabels_next_filtered) < 2:
                    print 'WARNING: not enought labels passed the size filter. increase template size or decrease size filter!'
                    continue

                distances = np.array(distances)
                nearestLabels_next_filtered = np.array(nearestLabels_next_filtered)
                indexes = np.argsort(distances)
                distances = distances[indexes].tolist()[:min(len(indexes),3)]
                nearestLabels_next_filtered = nearestLabels_next_filtered[indexes].tolist()[:len(distances)]

                if val == 1: # negative example
                    for i, label1 in enumerate(nearestLabels_next_filtered):
                        for label2 in nearestLabels_next_filtered[i+1:]:
                            negativeDivisions[t].append([[idx], [label1], [label2]])
                elif val == 2: # positive example
                    # assume that the closest 2 are children
                    positiveDivisions[t].append([[idx], [nearestLabels_next_filtered[0]], [nearestLabels_next_filtered[1]]])
                else:
                    raise Exception, "shouldn't get here"

    return positiveDivisions, negativeDivisions


def getUniqueFeatNames(features):
   result = []
   for v in features.values():
      for n in v:
         if n == 'None':
            continue
         if n not in result:
            result.append(n)

   return result

def getImage_at(fn, path, t, t_axis, ch_axis, ch):
   with h5py.File(fn, 'r') as f:   
      shape = f[path].shape
      slicing = [ slice(None), ] * len(shape)
      slicing[t_axis] = t
      slicing[ch_axis] = ch
      ds = f[path][tuple(slicing)]
   return ds

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
    


def getVigraFeatures(ds_raw, ds, features, ndim=None):
   r = ds_raw.squeeze()
   l = ds.squeeze()
   
   if ndim != None:
      if len(r.shape) != ndim:
         print 'WARNING: this image does not have the right direction, adding newaxis'
         for i in range(ndim-len(r.shape)):
            r = r[..., np.newaxis]
            l = l[..., np.newaxis]
   return vigra.analysis.extractRegionFeatures(r.astype(np.float32), l.squeeze().astype(np.uint32), features=features, ignoreLabel=0)
      
def getMaximalBoundingBox( bbPair, segmentsPair ):
   coordLength = len(bbPair[0]['Coord<Minimum >'][0])
   bb_min = [None, ] * coordLength
   bb_max = [None, ] * coordLength

   for idx, segments in enumerate(segmentsPair):
      for segment_id in segments:
         for dim in range(coordLength):            
            if bb_min[dim] is None or bb_min[dim] > bbPair[idx]['Coord<Minimum >'][segment_id][dim]:
               bb_min[dim] = bbPair[idx]['Coord<Minimum >'][segment_id][dim]
            if bb_max[dim] is None or bb_max[dim] < bbPair[idx]['Coord<Maximum >'][segment_id][dim]:
               bb_max[dim] = bbPair[idx]['Coord<Maximum >'][segment_id][dim]
   for i in range(len(bb_min)):
      if bb_max[i] - bb_min[i] < 2:
         bb_max[i] = bb_min[i] + 2
   return bb_min, bb_max

def getImage( img, boundingBox_min, boundingBox_max ):
   assert len(img.shape) == len(boundingBox_min)
   assert len(boundingBox_min) == len(boundingBox_max)
   slicing = []
   for idx in range(len(img.shape)):
      slicing.append(slice(boundingBox_min[idx], boundingBox_max[idx]+1))
   return img[slicing]

def relabel( volume, replace, defaultValue = 0 ):
   mp = np.arange(0,np.amax(volume)+1, dtype=volume.dtype)
   mp[1:] = defaultValue
   mp[replace.keys()] = replace.values()
   return mp[volume.astype(np.uint32)]

def relabelImage( img, positiveSegments, positiveSegments2=None, label=1, label2=2):
   replace = {}
   for v in positiveSegments:
      replace[v] = label
   if positiveSegments2 != None:
      for v in positiveSegments2:
         replace[v] = label2
   return relabel( img, replace )

def trainClassifier(featureArrays, labels):
   rf = vigra.learning.RandomForest(treeCount=100)
   featureArray = np.array(featureArrays).astype(np.float32)
   
   # replace nans
   rows, cols = np.where(np.isnan(featureArray) + np.isinf(featureArray))
   idx = (rows, cols)
   rows = list(set(rows.flat))
   cols = list(set(cols.flat))
   featureArray[idx] = 0

   labels = np.array(labels)[:,np.newaxis].astype(np.uint32)
   oob = rf.learnRF(featureArray, labels)
   return rf, oob
  
def writeToFile(out_fn, internal_path, rf, features, margin=None):
   rf.writeHDF5(out_fn, pathInFile=internal_path+'/rf')
   # write feature names to file
   with h5py.File(out_fn, 'r+') as f:
      if 'Features' in f[internal_path].keys():
         del f[internal_path]['Features']
      g = f[internal_path].create_group('Features')
      for key in getFeatureNamesOrdered(features.keys()):
         g.create_dataset(name=str(key), data=getFeatureNamesOrdered(features[key]))
      if margin != None:
         if 'margin' in f.keys():
            del f['margin']
         f.create_dataset(name='margin', data=margin)
      if 'timestamp' in f.keys():
         del f['timestamp']
      f.create_dataset(name='timestamp', data=str(datetime.date.today()))

def initFeatExtractors( features ):
   c = pgmlink.FeatureExtractorCollection()

   for operator in getFeatureNamesOrdered(features.keys()):
      for name in getFeatureNamesOrdered(features[operator]):
         c.addToVector(operator, name)

   return c

def makeTraxel( features, idx ):
   trax = pgmlink.Traxel()
   for name in features.keys():
      if 'Global' in name:
         continue
      val = features[name][idx]
      if isinstance(val, np.ndarray):
         val = val.flatten()
         trax.add_feature_array(name, len(val))
         for i in range(len(val)):
            trax.set_feature_value(name, i, float(val[i]))
      else: 
         trax.add_feature_array(name, 1)
         trax.set_feature_value(name, 0, float(val))
 
   return trax

def getFeatureNamesOrdered( feat_list ):
   return sorted(feat_list)

def getRegionFeatures(feat_extractor, feats_cur, idx_cur=1):   
   trax_cur = makeTraxel(feats_cur, idx_cur)
   
   return pgmlink.extractFeatures(feat_extractor, trax_cur).squeeze()

def getPairFeatures(feat_extractor, feats_cur, feats_next, idx_cur=1, idx_next=1):   
   trax_cur = makeTraxel(feats_cur, idx_cur)
   trax_next = makeTraxel(feats_next, idx_next)
   
   return pgmlink.extractFeatures(feat_extractor, trax_cur, trax_next).squeeze()

def getTripletFeatures(feat_extractor, feats_cur, feats_next, idx_cur=1, idx_next1=1, idx_next2=2):
   print len(feats_cur['Count']), idx_cur
   print len(feats_next['Count']), idx_next1, idx_next2
   trax_parent = makeTraxel(feats_cur, idx_cur)
   trax_child1 = makeTraxel(feats_next, idx_next1)
   trax_child2 = makeTraxel(feats_next, idx_next2)
   
   return pgmlink.extractFeatures(feat_extractor, trax_parent, trax_child1, trax_child2)


def getTransitionClassifier(raw_fn, pathRaw, fn, annotationPath, features, segmentImg_fn, segmentImg_path, t_axis, ch_axis, ch, out_fn, out_path, margin=None):
   print 'Transition Classifier:'
   print '======================'
   print '  reading annotations...'
   segment2Label, label2Segment = getAnnotations(fn, annotationPath)
   
   vigraFeatNames = getUniqueFeatNames(features)
   if 'None' in vigraFeatNames:
      vigraFeatNames.remove('None')
   if 'Count' not in vigraFeatNames:
      vigraFeatNames.append('Count')
   print '  initializing feature extractors...'
   featExtractor = initFeatExtractors(features).getExtractors()
   
   segmentBoundingBoxes = {}
   
   featureArrays = []
   # same ordering as featureArrays
   # label 1 means negative example, label 2 means positive example
   labels = []

   for t in label2Segment.keys():
      print '  processing timestep', t
      if t+1 not in label2Segment.keys():
         continue

      positivePairs_at = []
      negativePairs_at = []
      
      for label in label2Segment[t]:
         # if this label is not present in the next time step, this pair has already been added or is invalid
         if label not in label2Segment[t+1].keys():
            continue
         pair = [ label2Segment[t][label], label2Segment[t+1][label] ]
         if label < 0:
            negativePairs_at.append(pair)
         elif label > 0:
            positivePairs_at.append(pair)
         else:
            raiseException, "labels must not be zero"
      
      if len(positivePairs_at) + len(negativePairs_at) == 0:
         # nothing to do in this time step
         continue
      
      print '  reading segment images...'
      segmentImg_cur = getImage_at(segmentImg_fn, segmentImg_path, t, t_axis, ch_axis, ch).squeeze()
      segmentImg_next = getImage_at(segmentImg_fn, segmentImg_path, t+1, t_axis, ch_axis, ch).squeeze()
      rawImg_cur = getImage_at(raw_fn, pathRaw, t, t_axis, ch_axis, ch).squeeze()
      rawImg_next = getImage_at(raw_fn, pathRaw, t+1, t_axis, ch_axis, ch).squeeze()
      ndim = len(segmentImg_cur.squeeze().shape)

      #print '  computing bounding boxes...'
      #if t not in segmentBoundingBoxes.keys():      
      #   segmentBoundingBoxes[t] = getVigraFeatures(rawImg_cur, segmentImg_cur, ['Coord<Minimum >', 'Coord<Maximum >'])
      #segmentBoundingBoxes[t+1] = getVigraFeatures(rawImg_next, segmentImg_next, ['Coord<Minimum >', 'Coord<Maximum >'])
      regionIdxs = []
      relabelDict_cur = {}
      relabelDict_next = {}
      count = 1
      for idx, pair in enumerate(positivePairs_at + negativePairs_at):
         found = False
         for r in pair[0]:
            if r in relabelDict_cur.keys():
                found = True
                break
         for r in pair[1]:
            if r in relabelDict_next.keys():
                found = True
                break
         if found:
            print("Warning: Ignoring transition {} because regions have been used before".format(pair))
            regionIdxs.append(-1)
            continue
         regionIdxs.append(count)
         for r in pair[0]:
            relabelDict_cur[r] = count
         for r in pair[1]:
            relabelDict_next[r] = count
         count += 1
      seg_cur = relabel(segmentImg_cur, relabelDict_cur)
      seg_next = relabel(segmentImg_next, relabelDict_next)

      print '  extracting pairwise features...'
      feats_cur = getVigraFeaturesFromIlastik(rawImg_cur, seg_cur, vigraFeatNames, ndim, margin=margin)
      feats_next = getVigraFeaturesFromIlastik(rawImg_next, seg_next, vigraFeatNames, ndim, margin=margin)

      intersects = pgmlink.IntersectsDictionary()
      intersects.getIntersects(seg_cur.astype(np.uint32), seg_next.astype(np.uint32))

      for idx,pair in enumerate(positivePairs_at + negativePairs_at):
         #bb_min, bb_max = getMaximalBoundingBox( [segmentBoundingBoxes[t], segmentBoundingBoxes[t+1] ], pair )
         #seg_cur = getImage( segmentImg_cur, bb_min, bb_max )
         #seg_next = getImage( segmentImg_next, bb_min, bb_max )
         #img_cur = getImage( rawImg_cur, bb_min, bb_max )
         #img_next = getImage( rawImg_next, bb_min, bb_max )

         #seg_cur = relabelImage(seg_cur, pair[0])
         #seg_next = relabelImage(seg_next, pair[1])
         
         #feats_cur = getVigraFeatures(img_cur, seg_cur, vigraFeatNames, ndim)
         #feats_next = getVigraFeatures(img_next, seg_next, vigraFeatNames, ndim)

         if regionIdxs[idx] == -1:
             continue
         trax1 = makeTraxel(feats_cur, regionIdxs[idx])
         trax2 = makeTraxel(feats_next, regionIdxs[idx])

         intersects.addToTraxel(trax1, trax1.Id)
         
         if idx < len(positivePairs_at):
            labels.append(2)
         else:
            labels.append(1)
         featureArrays.append(pgmlink.extractFeatures(featExtractor, trax1, trax2).squeeze())
         #featureArrays.append(getPairFeatures(featExtractor, feats_cur, feats_next, regionIdxs[idx], regionIdxs[idx]))
   
   print '  training random forest...'      
   print "From labels:"
   print labels
   rf, oob = trainClassifier(featureArrays, labels)
   print '  =================================='
   print '  Transition Classifier: oob =', oob
   print '  =================================='

   print '  writing to file...'
   writeToFile(out_fn, out_path, rf, features, margin)
 

def getDivisionClassifier(raw_fn, pathRaw, fn, annotationPath, features, segmentImg_fn, segmentImg_path, t_axis, ch_axis, ch, out_fn, out_path, divisionPath, margin=None, from_ilastik=False):
   print 'Division Classifier:'
   print '===================='
   print '  initializing feature extractors...'
   featExtractor = initFeatExtractors(features).getExtractors()
   vigraFeatNames = getUniqueFeatNames(features)
   featureArrays = []
   labels = []
   if not from_ilastik:
       print '  reading annotations...'   
       segment2Label, label2Segment = getAnnotations(fn, annotationPath)
       positiveDivisions, negativeDivisions = getDivisionAnnotations(fn, divisionPath)
       
       if 'None' in vigraFeatNames:
           vigraFeatNames.remove('None')
       
       
       # same ordering as featureArrays
       # label 1 means negative example, label 2 means positive example

       positiveTriplets = {}
       negativeTriplets = {}
       timesteps = set()

       for div_idx, div in enumerate(positiveDivisions + negativeDivisions):
          assert len(div) == 3

          tParent = None
          for t in label2Segment.keys():
             if div[0] in label2Segment[t]:
                tParent = t
                break
          if tParent == None:
              print 'WARNING: parent is not contained in label set, continue...'
              continue
          assert tParent != None, 'the parent must be contained in the label set'
          assert tParent+1 in label2Segment.keys(), 'the label set must have another time step'
          assert div[1] in label2Segment[tParent+1].keys(), 'child1 must be in the label2Segment[t+1] set'
          assert div[2] in label2Segment[tParent+1].keys(), 'child2 must be in the label2Segment[t+1] set'
          
          triplet = [ label2Segment[tParent][div[0]], label2Segment[tParent+1][div[1]], label2Segment[tParent+1][div[2]] ]

          if tParent not in positiveTriplets:
             positiveTriplets[tParent] = []
          if tParent not in negativeTriplets:
             negativeTriplets[tParent] = []

          if div_idx < len(positiveDivisions):
             positiveTriplets[tParent].append(triplet)
          else:
             negativeTriplets[tParent].append(triplet)
          timesteps.add(tParent)
           
   else:
      print '  reading annotations...'   
      positiveTriplets, negativeTriplets = getDivisionTripletsIlastik(fn)
      timesteps = []
      #timesteps = set(positiveTriplets.keys()) | set(negativeTriplets.keys())
      for t in set(positiveTriplets.keys()) | set(negativeTriplets.keys()):
          if t in positiveTriplets.keys() and len(positiveTriplets[t]) > 0:
              timesteps.append(t)
              continue
          if t in negativeTriplets.keys() and len(negativeTriplets[t]) > 0:
              timesteps.append(t)

       #segmentBoundingBoxes = {}
   for t in timesteps:
      print '  timestep =', t
      print '  reading segment images...'
      segmentImg_cur = getImage_at(segmentImg_fn, segmentImg_path, t, t_axis, ch_axis, ch).squeeze()
      segmentImg_next = getImage_at(segmentImg_fn, segmentImg_path, t+1, t_axis, ch_axis, ch).squeeze()
      rawImg_cur = getImage_at(raw_fn, pathRaw, t, t_axis, ch_axis, ch).squeeze()
      rawImg_next = getImage_at(raw_fn, pathRaw, t+1, t_axis, ch_axis, ch).squeeze()
      ndim = len(segmentImg_cur.squeeze().shape)

      if not from_ilastik:
          #print '  computing bounding boxes...'
          #if t not in segmentBoundingBoxes.keys():      
          #   segmentBoundingBoxes[t] = getVigraFeatures(rawImg_cur, segmentImg_cur, ['Coord<Minimum >', 'Coord<Maximum >'])
          #segmentBoundingBoxes[t+1] = getVigraFeatures(rawImg_next, segmentImg_next, ['Coord<Minimum >', 'Coord<Maximum >'])
          regionIdxs = []
          relabelDict_cur = {}
          relabelDict_next = {}
          count = 1
          count_next = 1
          for idx, triplet in enumerate(positiveTriplets[t] + negativeTriplets[t]):
             found = False
             for r in triplet[0]:
                if r in relabelDict_cur.keys():
                    found = True
                    break
             for r in triplet[1]+triplet[2]:
                if r in relabelDict_next.keys():
                    found = True
                    break
             if found:
                print("Warning: Ignoring division {} because regions have been used before".format(triplet))
                regionIdxs.append([-1])
                continue
             for r in triplet[0]:
                relabelDict_cur[r] = count
             for r in triplet[1]:
                relabelDict_next[r] = count_next
             count_next += 1
             for r in triplet[2]:
                relabelDict_next[r] = count_next
             regionIdxs.append([count,count_next-1,count_next])
             count_next += 1            
             count += 1
          seg_cur = relabel(segmentImg_cur, relabelDict_cur)
          seg_next = relabel(segmentImg_next, relabelDict_next)
      else:
          seg_cur = segmentImg_cur
          seg_next = segmentImg_next
          regionIdxs = {}
          for idx, triplet in enumerate(positiveTriplets[t] + negativeTriplets[t]):
             regionIdxs[idx] = [triplet[0][0], triplet[1][0], triplet[2][0]]

      print '  extracting triplet features...'   
      feats_cur = getVigraFeaturesFromIlastik(rawImg_cur, seg_cur, vigraFeatNames, ndim, margin=margin)
      feats_next = getVigraFeaturesFromIlastik(rawImg_next, seg_next, vigraFeatNames, ndim, margin=margin)

      for idx, triplet in enumerate(positiveTriplets[t] + negativeTriplets[t]):
         #bb_min, bb_max = getMaximalBoundingBox( [segmentBoundingBoxes[t], segmentBoundingBoxes[t+1], segmentBoundingBoxes[t+1] ], triplet )
         #seg_cur = getImage( segmentImg_cur, bb_min, bb_max )
         #seg_next = getImage( segmentImg_next, bb_min, bb_max )
         #raw_cur = getImage( rawImg_cur, bb_min, bb_max )
         #raw_next = getImage( rawImg_next, bb_min, bb_max )

         #seg_cur = relabelImage(seg_cur, triplet[0])
         #seg_next = relabelImage(seg_next, triplet[1], positiveSegments2=triplet[2])
         
         #feats_cur = getVigraFeatures(raw_cur, seg_cur, vigraFeatNames, ndim)
         #feats_next = getVigraFeatures(raw_next, seg_next, vigraFeatNames, ndim)
         
         if regionIdxs[idx][0]==-1:
            continue

         if idx < len(positiveTriplets[t]):
            labels.append(2)
         else:
            labels.append(1)
    
         featureArrays.append(getTripletFeatures(featExtractor, feats_cur, feats_next, regionIdxs[idx][0], regionIdxs[idx][1], regionIdxs[idx][2]).squeeze())

   print '  training random forest...'      
   rf, oob = trainClassifier(featureArrays, labels)
   print '  =================================='
   print '  Division Classifier: oob =', oob
   print '  =================================='

   print '  writing to file...'
   writeToFile(out_fn, out_path, rf, features, margin)


def getRegionClassifier(raw_fn, pathRaw, fn, annotationPath, features, segmentImg_fn, segmentImg_path, t_axis, ch_axis, ch, out_fn, out_path, margin=(2,2,2)):
   print 'Region Classifier'
   print '================='
   print '  reading annotations...'
   segment2Label, label2Segment = getAnnotations(fn, annotationPath)
   
   vigraFeatNames = getUniqueFeatNames(features)
   print '  initializing feature extractors...'
   featExtractor = initFeatExtractors(features).getExtractors()
   
   segmentBoundingBoxes = {}
   
   featureArrays = []
   # same ordering as featureArrays
   # label 1 means negative example, label 2 means positive example
   labels = []

   for t in label2Segment.keys():
      print '  processing timestep', t

      positiveRegions_at = []
      negativeRegions_at = []
      
      for label, region in label2Segment[t].items():
         if label < 0:
            negativeRegions_at.append(region)
         elif label > 0:
            positiveRegions_at.append(region)
         else:
            raiseException, "labels must not be zero"
      
      if len(positiveRegions_at) + len(negativeRegions_at) == 0:
         # nothing to do in this time step
         continue
      
      print '  reading segment image...'
      segmentImg_cur = getImage_at(segmentImg_fn, segmentImg_path, t, t_axis, ch_axis, ch).squeeze()
      rawImg_cur = getImage_at(raw_fn, pathRaw, t, t_axis, ch_axis, ch).squeeze()
      ndim = len(segmentImg_cur.squeeze().shape)

      #print '  computing bounding boxes...'
      #if t not in segmentBoundingBoxes.keys():      
      #   segmentBoundingBoxes[t] = getVigraFeatures(rawImg_cur, segmentImg_cur, ['Coord<Minimum >', 'Coord<Maximum >'])
      regionIdxs = []
      relabelDict = {}
      count = 1
      for idx, region in enumerate(positiveRegions_at + negativeRegions_at):
         found = False
         for r in region:
            if r in relabelDict.keys():
                found = True
                break
         if found:
            print("Warning: Ignoring region {} because regions have been used before".format(region))
            regionIdxs.append(-1)
            continue
         regionIdxs.append(count)
         for r in region:
            relabelDict[r] = count
         count += 1
      seg_cur = relabel(segmentImg_cur, relabelDict)      

      feats_cur = getVigraFeaturesFromIlastik(rawImg_cur, seg_cur, vigraFeatNames, ndim, margin=margin)
      
      print '  extracting region features...'
      for idx, region in enumerate(positiveRegions_at + negativeRegions_at):

         #bb_min, bb_max = getMaximalBoundingBox( [ segmentBoundingBoxes[t] ], [ region ] )
         #seg_cur = getImage( segmentImg_cur, bb_min, bb_max )
         #raw_cur = getImage( rawImg_cur, bb_min, bb_max )

         #seg_cur = relabelImage(seg_cur, region)
         
         #feats_cur = getVigraFeaturesFromIlastik(raw_cur, seg_cur, vigraFeatNames, ndim, margin=margin)
         
         if regionIdxs[idx] == -1:
            continue

         if idx < len(positiveRegions_at):
            labels.append(2)
         else:
            labels.append(1)

         featureArrays.append(getRegionFeatures(featExtractor, feats_cur,idx_cur=regionIdxs[idx]))
         #assert len(featureArrays[0]) == len(featureArrays[-1])
   
   print '  training random forest...'      
   rf, oob = trainClassifier(featureArrays, labels)
   print '  =================================='
   print '  Region Classifier: oob =', oob
   print '  =================================='

   print '  writing to file...'
   writeToFile(out_fn, out_path, rf, features, margin)



if __name__ == '__main__':
   import sys
   argv = sys.argv

   parser = optparse.OptionParser(usage="%prog [options]")  
   parser.add_option('--output', '-o', dest='out_fn', default=None, type=str, help='classifier output file (hdf5) [default=%default]')

   parser.add_option('--oversegmentation', dest='segmentImg_fn', default=None, type=str, help='input file (hdf5) for oversegmentation [default=%default]')
   parser.add_option('--internalOversegmentation', dest='segmentImg_path', default='exported_data', type=str, help='internal hdf5 path to oversegmentation dataset [default=%default]')
   
   parser.add_option('--raw', '-r', dest='raw_fn', default=None, type=str, help='h5 file of raw data; keep None if watershed should be computed on prediction map rather than on the raw data [default=%default]')
   parser.add_option('--internalPathRaw', dest='pathRaw', default='exported_data', type=str, help='internal hdf5 path to raw data [default=%default]')

   parser.add_option('--labelProject', '-p', dest='fn', default=None, type=str, help='ilastik project containing the labels [default=%default]')
   parser.add_option('--internalPathTransitionLabels', dest='pathTransitions', default='/ManualTransitions/Labels/0', type=str, help='internal path of transition labels [default=%default]')
   parser.add_option('--internalPathRegionLabels', dest='pathRegions', default='/ManualRegions/Labels/0', type=str, help='internal path of region labels [default=%default]')
   parser.add_option('--internalPathDivisionLabels', dest='pathDivisions', default='/ManualDivisions/Labels/0', type=str, help='internal path of division labels [default=%default]')
   
   parser.add_option('--t-axis', dest='t_axis', default=0, type=int, help='index of time-axis [default=%default]')
   parser.add_option('--ch-axis', dest='ch_axis', default=-1, type=int, help='index of channel-axis [default=%default]')
   parser.add_option('--segmentation-channel', dest='ch', default=0, type=int, help='channel of segmentation [default=%default]')

   options, args = parser.parse_args()
   if options.fn == None:
      parser.print_help()
      sys.exit(0)

   div_from_ilastik = False

   neighbor_feats = ['Mean', 'Variance', 'Skewness', 'Kurtosis', 'Maximum', 'Minimum']
   neighbor_feats_neigh = [ x + str(',') + x + str(' in neighborhood') for x in neighbor_feats ]
   neighbor_feats_obj = [ x + str(',') + x + str(' in object and neighborhood') for x in neighbor_feats ]

   transitionFeatures = { 'AbsDiff': ['RegionCenter', 'Count', 'Mean', 'Variance', 'Sum', 'Coord<ArgMaxWeight>', 
                                      'Coord<ArgMinWeight>', 'Weighted<RegionCenter>', 'RegionRadii',
                                      'Skewness', 
                                      'Quantiles', 'Kurtosis', 'Maximum', 'Minimum',
                                      'RegionAxes', 'Weighted<RegionRadii>', 'Weighted<RegionAxes>',
                                    ],
                          'Ratio': ['Count', 'Mean', 'Sum',  
                                    'Variance', 'Skewness', 'Kurtosis', 'Maximum', 'Minimum',
                                    'RegionRadii', 'Quantiles',  
                                    'RegionAxes', 'Weighted<RegionRadii>', 'Weighted<RegionAxes>'
                                 ],
                           '_Ratio': ['Minimum,Maximum','Maximum,Minimum'] +
                                     neighbor_feats_obj +
                                     neighbor_feats_neigh,
                          'IntersectionUnionRatio' : ['None']
                       }
   # transitionFeatures = { 'AbsDiff': ['RegionCenter', 'Count', 'CenterOfMass'],
   #                        'Ratio': ['Count',],
   #                        'IntersectionUnionRatio' : ['Count']
   #                      }
   
   divisionFeatures = { #'ParentRatio': ['Mean', 'Count', 'Variance', 'Sum'],
                        #'ParentSquaredDifferencesRatio': ['RegionCenter', 'Mean', 'Count', 'Variance', 'Sum'],
                        'MaxParentRatio': transitionFeatures['Ratio'],
                        'MinParentRatio': transitionFeatures['Ratio'],
                        'MeanParentRatio': transitionFeatures['Ratio'],
                        'MinParentSquaredDifference': transitionFeatures['AbsDiff'],
                        'MaxParentSquaredDifference': transitionFeatures['AbsDiff'],
                        'MeanParentSquaredDifference': transitionFeatures['AbsDiff'],
                        'RatioParentSquaredDifference': transitionFeatures['AbsDiff'],
                       }
                     
   regionFeatures = { 'Identity': ['RegionRadii', 'Mean', 'Count', 'Variance', 'Sum', 'Skewness',
                                    'Quantiles', 'Kurtosis', 'Maximum', 'Minimum',
                                    'RegionAxes', 'Weighted<RegionRadii>', 'Weighted<RegionAxes>',
                                    ] + [x + str(" in neighborhood") for x in neighbor_feats] +
                                        [x + str(" in object and neighborhood") for x in neighbor_feats],
                       '_Ratio': neighbor_feats_obj + neighbor_feats_neigh
                     }
                      
   neighborhood_margin = (2,2,2)

   getTransitionClassifier(options.raw_fn, options.pathRaw, options.fn, options.pathTransitions, transitionFeatures, options.segmentImg_fn, options.segmentImg_path, options.t_axis, options.ch_axis, options.ch, options.out_fn, 'Transitions')
   getDivisionClassifier(options.raw_fn, options.pathRaw, options.fn, options.pathDivisions, divisionFeatures, options.segmentImg_fn, options.segmentImg_path, options.t_axis, options.ch_axis, options.ch, options.out_fn, 'Divisions', '/ManualDivisions/Divisions/0', from_ilastik=div_from_ilastik)
   getRegionClassifier(options.raw_fn, options.pathRaw, options.fn, options.pathRegions, regionFeatures, options.segmentImg_fn, options.segmentImg_path, options.t_axis, options.ch_axis, options.ch, options.out_fn, 'Regions', margin=neighborhood_margin)

