#ifndef FILE_READER_H
#define FILE_READER_H

#include <string>
#include <vector>
#include <map>
#include <boost/shared_array.hpp>

namespace Tracking {

// General interface class for reading a data file with the features 
// required for the tracking
//
class FileReader {
public:

// Construct from file name and Z / XY resolution ratio
  FileReader( const std::string& fileName, double xyzRatio=1. );

// Read the input features from the encapsulated data file. 
// Upon successful completion, nCells contains the number of cells,
// nFeats contains the number of features and feats contains the nFeats*nCells
// feature data (first all features for cell no. 1 stacked behind each other,
// then all features for cell no. 2, ...).
  int readInputs(int& nCells, int& nFeats, 
		 boost::shared_array<double>& feats);
  
// Read-only access to file name
  std::string getFileName() const;

// Get a vector of all valid cell labels for this file
  void getValidLabels(std::vector<int>& validLabels) const;

// Get the Z / XY resolution ratio
  double getXYZRatio() const;
protected:
  virtual int readInputsSpecific(int& nCells, int& nFeats,
                                 boost::shared_array<double>& feats) = 0;
// Default implementation: returns all numbers between 1 and nCells
  virtual void getValidLabelsSpecific(std::vector<int>& validLabels) const;
private:
  std::string fileName_;
  int nCells_;
// Z to XY resolution ratio (should be >= 1.)
  double xyzRatio_; 
};

// Read the files from the EMBL Philipp Keller datafiles: see
// http://www.embl.de/digitalembryo/

class FileReaderEmbl : public FileReader {
public:
// Constructor
  FileReaderEmbl( const std::string& fileName, double xyzRatio=1. );
protected:   
// Read the input coordinates from a specified HDF5 file in the 
// EMBL format (see http://www.embl.de/digitalembryo/)
   virtual int readInputsSpecific(int& nCells, int& nFeats,
                                 boost::shared_array<double>& feats);
};

// Reader for the input coordinates from the HCI format for digital embryo
// files (see 
// https://gorgonzola.iwr.uni-heidelberg.de/intern/wiki/index.php/Interface )
class FileReaderHci : public FileReader {
public:
// Constructor with initialization of admissible named features
  FileReaderHci( const std::string& fileName, double xyzRatio=1. );
// Determine which features shall be used for the tracking
  int setFlag( unsigned long flag );
protected:
  virtual int readInputsSpecific( int& nCells, int& nFeats, 
                                  boost::shared_array<double>& feats );
  virtual void getValidLabelsSpecific( std::vector<int>& validLabels ) const;
private:
  // Named Feature structure that contains all necessary information for
  // retrieving a feature and the corresponding information from an HDF5 file
  struct NamedFeature {
    std::string name_; // e.g. "position"
    int nSubFeats_;    // how many subfeatures this feature comprises, e.g.
                       // 3 for a position
    boost::shared_array<double> weights_; // weights to be used for this feature
  };
  // Named flags that can be used to decide which data to extract from
  // the file
  // key (long) -- Flag used in the "featureconfig" data element and showing
  //               which features have been computed
  // value (NamedFeature) -- Corresponding feature to be stored             
  std::map<unsigned long, NamedFeature> namedFlags_;
    
  unsigned long flag_; // flag deciding which features shall be extracted

  // Vector containing 1's for valid cell labels and 0's for invalid labels
  std::vector<char> labelValidity_;
};


}

#endif /* FILE_READER_H */
