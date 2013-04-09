#!/usr/bin/python

import vigra
from pgmlink import feature_array, gmm_priors_and_centers
import numpy as np
import sys

def append_np_array_to_feature_array(np_array, feat_array):
    for el in np_array:
        feat_array.push_back(el)

if __name__ == "__main__":

    for weight in range(1,11):
        w = 1.0/weight
        print "weight =", w
        data = feature_array()
        # cluster around 0
        c = np.random.randn(50)
        append_np_array_to_feature_array(c, data)
        # cluster around off1
        off1=6
        c = np.random.randn(50) + off1
        append_np_array_to_feature_array(c, data)
        # cluster around off2
        off2=18
        c = np.random.randn(50) + off2
        append_np_array_to_feature_array(c, data)
        centers = feature_array()
        priors = feature_array()
        n = 1
        k_max = 5
        gmm_priors_and_centers(data, priors, centers, k_max, n, w)
        for k in range(k_max):
            print "   -> k = %d:" % (k+1,), priors[k]
        print
    print

    for weight in range(1,11):
        w = 1.0/weight
        print "weight =", w
        # cluster around 0
        c = np.random.randn(50)
        # cluster around off1
        off1=6
        c = np.append(c,np.random.randn(50) + off1)
        # cluster around off2
        off2=18
        c = np.append(c,np.random.randn(50) + off2)
        centers = feature_array()
        priors = feature_array()
        n = 1
        k_max = 5
        c = c.reshape((1,) + c.shape)
        gmm_priors_and_centers(c, priors, centers, k_max, n, w)
        for k in range(k_max):
            print "   -> k = %d:" % (k+1,), priors[k]
        print centers[0], centers[1], centers[2], centers[3]
        print
