import numpy as np
import vigra
import pgmlink


def relabel( volume, replace ):
    mp = np.arange(0,np.amax(volume)+1, dtype=volume.dtype)
    mp[1:] = 0
    mp[replace.keys()] = replace.values()
    return mp[volume]


# goes through segment_image and builds up a region adjacency graph
# where neighboring segments are connected by a weighted edge 
# returns this graph
def regionAdjacencyGraphFromSegmentImage(segmentImage):
   graph = pgmlink.RegionAdjacencyGraph()
   graph.buildGraph(segmentImage.squeeze().astype(np.float))
   assert pgmlink.countNodes(graph) == np.max(segmentImage)
   return graph

# copies the provided graph, then
# merges all nodes in the given graph where the edge nodes are above
# the provided threshold, then returns this graph
def mergeNodesThreshold(graph, percentile=0.8, greedy=True):
   nnodes = pgmlink.countNodes(graph)
   nedges = pgmlink.countEdges(graph)
   print '    number of nodes in graph:', nnodes
   print '    number of edges in graph:', nedges
   weights = list(graph.getEdgeWeightMap().values())

   if len(weights) > 0:
       threshold = np.percentile(weights, percentile*100)
       print '    threshold =', threshold
       graph.mergeNodesThreshold(threshold, greedy)

       print '    contracted', nedges-pgmlink.countEdges(graph), 'edges and removed/contracted', \
           nnodes-pgmlink.countNodes(graph), 'nodes'
       print
   return graph

# returns a labelimage of the given shape where each node in 
# the region graph has its unique identifier
def labelImageFromRegionGraph(graph, segmentImage, alreadyAddedRegions=[]):
   labelsVector = graph.getLabelsVector()
   regions = graph.getRegions()
   reg_ids = []
   reg_labels = []
   for r in regions:
      reg_ids.append(r.id)
      reg_labels.append(list(r.contains_labels))

   replace = {}
   for idx, labels in enumerate(labelsVector):
      labels = sorted(list(labels))
      reg_id = None
      for i, v in enumerate(reg_labels):
         if labels == v:
            reg_id = reg_ids[i]
            break
      if reg_id in alreadyAddedRegions:
         for l in labels:
            replace[l] = 0
         continue
      alreadyAddedRegions.append(reg_id)

      for l in labels:
         if reg_id is None:
            raise Exception, "could not find reg_id"
         else:
            replace[l] = reg_id

   result = relabel(segmentImage, replace)

   return result

def getConflictSets(graph):
   conflictsMap = graph.getConflictSets()
   return conflictsMap

def getConnectedComponents(graph):
   cc = graph.getConnectedComponentIds()
   return cc

def getRegionToConnectedCompMap(graph):
   regions = graph.getRegions()
   cc = list(getConnectedComponents(graph))

   regionToConnectedComponent = {}

   cc_labels = []
   cc_ids = []   


   for r in regions:
      if r.id in cc:
         assert int(r.id) not in regionToConnectedComponent.keys()
         regionToConnectedComponent[int(r.id)] = int(r.id)
         cc_ids.append(int(r.id))
         cc_labels.append(list(r.contains_labels))
   
   for r in regions:
      if int(r.id) in regionToConnectedComponent.keys():
         continue
      found = False
      r_label = list(r.contains_labels)[0]
      for i, cc_l in enumerate(cc_labels):
         if r_label in cc_l:
            regionToConnectedComponent[int(r.id)] = cc_ids[i]
            found = True
            break
      assert found, "the region must be contained in some connected component!"

   assert len(regionToConnectedComponent) == len(regions)

   return regionToConnectedComponent
   
def getConnectedComponentImage(regionMaps, segmentImage):
   replace = {}
   max_seg_id = np.max(segmentImage)
   for seg_id in range(1,max_seg_id):
      replace[seg_id] = regionMaps[seg_id]
   return relabel(segmentImage, replace)
   

def readImage(fn):
   return vigra.impex.readImage(fn).astype(np.uint16)


def writeImage(img, fn):
   vigra.impex.writeImage(img, fn)



if __name__ == "__main__":
   levels = 6

   print 'reading in segment image...'
   segmentImage = readImage("./segmentImage.tif")
   
   labelImages = []
   print 'building region adjacency graph...'   
   graph = regionAdjacencyGraphFromSegmentImage(segmentImage)
   labelImages.append(labelImageFromRegionGraph(graph, segmentImage))

   for level in range(levels-1):
      print 'merging regions for level', level, '...'
      graph = mergeNodesThreshold(graph, percentile=0.8-0.1*level) # 0.8, 0.7, 0.6,...
      print 'drawing labelImage for level', level, '...'
      labelImages.append(labelImageFromRegionGraph(graph, segmentImage))
   print 'get remaining connected components...'
   graph = mergeNodesThreshold(graph, percentile=0)
   print 'drawing labelImage for connected components...'
   labelImages.append(labelImageFromRegionGraph(graph, segmentImage))

   print 'writing label images...'
   for level in range(len(labelImages)):
      writeImage(labelImages[level], './labelImage_level' + str(level) + '.tif')

   print 'done'

   
   cc = getConnectedComponents(graph)
   cs = getConflictSets(graph)

   for c in cc:
      #print 'c =', c
      for vec in cs[c]:
         s = ''
         for v in vec:
            s += str(v) + ', '
         #print s

