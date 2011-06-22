#include "hdf5.h"
#include "hdf5_hl.h"
#include "ANN.h"
#include <boost/shared_array.hpp>
#include "lp_lib.h"
#include <iostream>

int main(){
  lprec *lp = make_lp(10,0);
  hid_t fid = H5Fcreate( "test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  boost::shared_array<int> data( new int[4] );
  hsize_t dims[2] = {2, 2};
  H5LTmake_dataset_int( fid, "/Test", 2, dims, data.get() );
  H5Fclose( fid );
  unsigned hdf5Vers[3];
  H5get_libversion(hdf5Vers,hdf5Vers+1,hdf5Vers+2);
  std::cout << "HDF5 version = " << hdf5Vers[0] << "." << hdf5Vers[1] << "." << hdf5Vers[2] << std::endl;
  ANNpointArray pts = annAllocPts( 5 , 3 );
  annDeallocPts( pts );
  delete_lp( lp );
  std::cout << "Library test ended successfully" << std::endl;
  return 0;
}

