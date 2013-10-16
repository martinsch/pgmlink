# coding: utf-8
import pgmlink as pgm
import numpy as np
c = pgm.FeatureExtractorCollection()
c3 = pgm.FeatureExtractorCollection()

c.addToVector('Ratio', 'Count')
c.addToVector('Ratio', 'Intensity')
c.addToVector('SquaredDiff', 'RegionCenter')
c.addToVector('AbsDiff', 'RegionCenter')
c.addToVector('SqrtSquaredDiff', 'RegionCenter')
c.addToVector('AbsDiff', 'Intensity')
c.addToVector('SquaredDiff', 'Count')
c3.addToVector('Ratio', 'Intensity')
c3.addToVector('MaxParentRatio', 'Count')
c3.addToVector('MinParentRatio', 'Count')
c3.addToVector('MeanParentRatio', 'Count')
c3.addToVector('MaxParentRatio', 'Count')
c3.addToVector('MinParentRatio', 'Count')
c3.addToVector('MeanParentRatio', 'Count')
c3.addToVector('MaxParentSquaredDifference', 'RegionCenter')
c3.addToVector('MinParentSquaredDifference', 'RegionCenter')
c3.addToVector('MeanParentSquaredDifference', 'RegionCenter')
c3.addToVector('RatioParentSquaredDifference', 'RegionCenter')

t = [None] * 3
t[0] = pgm.Traxel()
t[1] = pgm.Traxel()
t[2] = pgm.Traxel()
com = np.arange(9, dtype=np.float64).reshape((3,3))
count = np.arange(10,13, dtype=np.float64).reshape((3,1))
intensity = 50 - np.arange(22,25, dtype=np.float64).reshape((3,1))

for idx, trax in enumerate(t):
    trax.add_feature_array('RegionCenter', 3)
    trax.add_feature_array('Count', 1)
    trax.add_feature_array('Intensity', 1)
    for i in range(3):
        trax.set_feature_value('RegionCenter', i, com[idx, i])
    trax.set_feature_value('Count', 0, count[idx,0])
    trax.set_feature_value('Intensity', 0, intensity[idx,0])

ex = c.getExtractors()
res = pgm.extractFeatures(ex, t[0], t[1])
print res

res = pgm.extractFeatures(ex, t[1], t[2])
print res

res = pgm.extractFeatures(ex, t[0], t[2])
print res

try:
    res = pgm.extractFeatures(ex, t[0], t[1], t[2])
    print res
except RuntimeError, e:
    print "  Fails es expected:"
    print " ", e

ex = c3.getExtractors()
res = pgm.extractFeatures(ex, t[0], t[1], t[2])
print res

# res = pgm.calculateFeature(c.getCalculators(), feats, feats, feats)


# print res
