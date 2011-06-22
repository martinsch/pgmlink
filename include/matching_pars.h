#ifndef MATCHING_PARS_H
#define MATCHING_PARS_H

// Parameters for the matching process 
class MatchingPars {
public:
  // number of nearest neighbors to be considered as children
  int knn; 
  // maximum distance between a cell and its potential children
  double maxDist; 
  // cost for a split
  double cSplit; 
  
  double cDeath; // cost for a cell death
  
  double cAppear; // cost for an appearance

  double xyzRatio; // Ratio of lateral resolution to axial resolution (typically > 1)
  
  unsigned long flagFeatures; // flag determining which features shall be used
  
  MatchingPars(); // Constructor with default arguments
  
  int readOptions(int nOpts, char* opts[]); // read parameters from string array
};

#endif /* MATCHING_PARS_H */

