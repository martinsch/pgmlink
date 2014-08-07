#!/usr/bin/env python

import vigra
import numpy as np
import sys
import random
import h5py

def do_oversegment(img, sigma=0.3, method='dtSeededGradientMagnitude', withOpening=False, convert=False, seedSigma=3.0, seedDiscRadius=2, sizeFilterFrom=2):
    is3D = False
    if len(img.shape) == 3 and img.shape[-1] != 1:
       is3D = True

    seeds = None
    bg_seed = None

    # delete all small objects and fill holes because a (dilated) seed on it would seed the background as well...
    if withOpening:
      binary = vigra.filters.discOpening(img.astype(np.bool).astype(np.uint8), 3)
    else:
      binary = img.astype(np.bool).astype(np.uint8)

    print '    generating seeds...'
    # make a distance transform to the background and look for local maxima
    if method=='gradientSeededGradientMagnitude':
      dt = vigra.filters.gaussianGradientMagnitude(img.astype(np.float32), seedSigma) # scale can be adjusted to get different strengths of oversegmentation
      if is3D:
         seeds = vigra.analysis.extendedLocalMinima3D(dt.astype(np.float32)) 
      else:
         seeds = vigra.analysis.extendedLocalMinima(dt.astype(np.float32)) 
    elif method=='dtSeededGradientMagnitude' or method=='dtSeededDt':
      binary_aug = np.zeros(np.array(binary.shape)+2)
      slicing = [slice(None),] * len(binary_aug.shape)
      for i in range(len(slicing)):
         if binary_aug.shape[i] != 1:
            slicing[i] = slice(1, -1)
      slicing = tuple(slicing)
      binary_aug[slicing] = binary
      if is3D:
         dt = vigra.filters.distanceTransform3D(binary_aug.astype(np.float32), background=False)
         dt = dt[slicing]
         seeds = vigra.analysis.extendedLocalMaxima3D(dt.astype(np.float32)) # neighborhood can be 4 or 8 (default)
      else:
         dt = vigra.filters.distanceTransform2D(binary_aug.astype(np.float32), background=False, norm=2) # norm can be 0 (infinity norm), 1 (L1 norm) or 2 (Euclidean norm)
         dt = dt[slicing]
         seeds = vigra.analysis.extendedLocalMaxima(dt.astype(np.float32)) # neighborhood can be 4 or 8 (default)
    else:
      raise Exception, 'method not implemented'

    # the seeds should be a bit bigger than just 1 pixel, so dilate with a disc (no dilation, if radius == 1)
    if is3D:
        seeds = vigra.filters.multiBinaryDilation(seeds.astype(np.uint8), seedDiscRadius)
    else:
        seeds = vigra.filters.discDilation(seeds.astype(np.uint8), seedDiscRadius) # radius can be everything!

    # set the seeds on the background to 0
    seeds *= binary.astype(np.bool)

    # now every connected component is one seed
    if is3D:
      seeds = vigra.analysis.labelVolumeWithBackground(seeds, neighborhood=6, background_value=0).astype(np.uint32)
      
    else:
      seeds = vigra.analysis.labelImageWithBackground(seeds, neighborhood=4, background_value=0).astype(np.uint32)

    # for the background, sample random 0-intensity-value pixels and set a background seed
    bg_seed = np.max(seeds) + 1
    zeros = np.where(dt == 0)
    random.seed(42) # seed the random generator for reproducible experiments
    for i in range(1000):
      idx = random.randint(0,len(zeros[0])-1)
      if is3D:
         seeds[zeros[0][idx],zeros[1][idx],zeros[2][idx]] = bg_seed
      else:
         seeds[zeros[0][idx],zeros[1][idx]] = bg_seed

    # make sure that every connected background component has at least one seed
    inverseBinary = 1 - binary.astype(np.bool)
    if is3D:
      labelImg_background = vigra.analysis.labelVolumeWithBackground(inverseBinary.astype(np.uint8)).astype(np.uint32)
    else:
      labelImg_background = vigra.analysis.labelImageWithBackground(inverseBinary.astype(np.uint8)).astype(np.uint32)
    cc_background = np.unique(labelImg_background).tolist()
    for cc in cc_background:
       cc_zeros = np.where(cc_background == cc)
       idx = random.randint(0,len(cc_zeros[0]))
       if is3D:
          seeds[zeros[0][idx],zeros[1][idx],zeros[2][idx]] = bg_seed
       else:  
          seeds[zeros[0][idx],zeros[1][idx]] = bg_seed


    # compute the graussian gradient magnitude image
    print '    computing Gaussian gradient magnitude...'
    if method != 'dtSeededDt':
        img = vigra.filters.gaussianGradientMagnitude(img.astype(np.float32), sigma)
    else:
        img = (np.max(dt) - dt)*binary
        
    # run a watershed, seeded if seeds are given, otherwise local maxima
    print '    computing watershed...'
    
    res = vigra.analysis.watershedsNew(img,seeds=seeds)
    clusters = np.unique(res[0]).tolist()
    print '    ' + str(len(clusters)) + ' clusters extracted.'
    if bg_seed is not None and bg_seed in clusters:
       clusters.remove(bg_seed)
    clusters = np.array(clusters)
    if convert:
       outChannels = 3
       dtype = np.uint8
    else:
       outChannels = 1
       count = 1
       dtype = np.uint16
    if is3D:
       out = np.zeros(img.shape[:3] + (outChannels,), dtype=dtype)
    else:
       out = np.zeros(img.shape[:2] + (outChannels,), dtype=dtype)
    for c in clusters:        
        ind = (res[0] == c)
        if ind.shape[-1] == 1:
            ind = ind[...,0]
        if sizeFilterFrom != None and np.count_nonzero(ind) <= sizeFilterFrom:
            continue
        for i in xrange(out.shape[-1]):
            if convert:
               out[...,i][ind] = np.random.randint(0,255)
            else:
               out[...,i][ind] = count
               count += 1
    return res, out


def readHdf5dataset(fn, internal_path, t, ch, t_axis=0, ch_axis=-1):
    with h5py.File(fn, 'r') as f:
        ds = f[internal_path]
        slicing = [slice(None),] * len(ds.shape)
        slicing[t_axis] = t
        slicing[ch_axis] = ch
        img = ds[tuple(slicing)].astype(np.float32)
    return img


def readTif2D(fn):
    img = vigra.impex.readImage(fn)
    return img

def writeTif2D(out, out_fn):
    vigra.impex.writeImage(out, out_fn)

def writeHdf5(out, out_fn, internal_path, original=None):
    with h5py.File(out_fn, 'w') as f:
        ds = f.create_dataset(internal_path, data=out, dtype=np.uint32)
        if original is not None:
           if 'original' in f.keys():
              del f['original']
           orig = f.create_dataset(name='original', data=original)

def binarize(img, threshold):
    return (img>=float(threshold)).astype(np.uint8)

def multiplyWith(img1, img2, dtype=np.uint8):
    return (img1 * img2).astype(dtype)

def normalize(img):
    maximum = np.max(img)
    a = 255./maximum
    res = (a * img).astype(np.uint8)
    assert np.max(res) == 255
    assert np.min(res) == 0
    return res
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("watershed")
    parser.add_argument('--predictionMap', '-p', required=True)
    parser.add_argument('--raw', '-r', default=None, type=str)
    parser.add_argument('--sigma', '-s', default=0.0, type=float)
    parser.add_argument('--input', '-i', default='dtSeededGradientMagnitude')
    parser.add_argument('--internalPathPredMap', default='exported_data')
    parser.add_argument('--internalPathRaw', default='/volume/data')
    parser.add_argument('--threshold', default=128, type=float)
    parser.add_argument('-t', default=0, type=int)
    parser.add_argument('-ch', default=1, type=int)
    parser.add_argument('--opening', '-o', action='store_true', default=False)
    parser.add_argument('--normalize', action='store_true', default=False)
    args = vars(parser.parse_args())
    fn = args['predictionMap']
    out_fn = fn
    if args['raw'] is not None:
       out_fn = args['raw']
    out_fn = ''.join(['.'.join(out_fn.split('.')[:-1]),] + ['_normalize=', str(args['normalize']), '_binarize-threshold=', str(args['threshold']), '_watershed=', args['input'], 'sigma=', str(args['sigma']), ',opening=', str(args['opening']), '.'] + [fn.split('.')[-1]])

    isH5 = False
    if fn.split('.')[-1] in ['h5','ilp','hdf5']:
       print 'WARNING: multiple time steps are not implemented yet'
       predMap = readHdf5dataset(fn, args['internalPathPredMap'], args['t'], args['ch'])
       isH5 = True
    else:
       predMap = readTif2D(fn)

    predMap = normalize(predMap)
    binaryPredMap = binarize(predMap, args['threshold'])
    if args['raw'] is None:
       img = predMap
    else:
       img = readHdf5dataset(args['raw'], args['internalPathRaw'], args['t'], 0)
       # if the axis order does not match, swap first and last axis
       if img.shape != predMap.shape:
          img = np.swapaxes(img,0,2)
          assert img.shape == predMap.shape
       if img.shape != predMap.shape:
          raise RuntimeError, "shapes do not match"
    
   
    img = multiplyWith(img, binaryPredMap, dtype=img.dtype)
    res, out = do_oversegment(img,
                         args['sigma'],
                         method='gradientSeededGradientMagnitude',
                         withOpening=args['opening'],
                         convert=True, # convert to RGB
                         )
    
    print out_fn
    if isH5:
       writeHdf5(out, out_fn, 'watershed_result', original=img[...,np.newaxis])
    else:
       writeTif2D(out, out_fn)
