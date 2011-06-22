#include "hdf5.h"
#include "hdf5_hl.h"
#include "ANN.h"
#include <boost/shared_array.hpp>
#include "lp_lib.h"
#include <ilcplex/cplex.h>

#include <iostream>

void reportCplexError(CPXENVptr env, int errcode){
  char errbuffer[4096];
  CPXCCHARptr outptr = CPXgeterrorstring(env, errcode, errbuffer);
  if( outptr!=NULL ){
    fprintf(stderr, "%s\n", errbuffer);
  } else {
    fprintf(stderr, "Unknown CPLEX error flag: %d\n", errcode); 
  }
}


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
  int status;
  CPXENVptr env = CPXopenCPLEX(&status);
  if( env==NULL ){
    reportCplexError( env, status );
    std::cerr << "Error while initializing CPLEX environment" << std::endl;
    return 1;
  } 
  CPXLPptr lpCpl = CPXcreateprob(env, &status, "Problem 1");
  if( lpCpl == NULL ){
    reportCplexError(env, status);
    std::cerr << "Error while creating CPLEX LP problem" << std::endl;
  } else {
    status = CPXfreeprob(env, &lpCpl );
    if( status ){
      std::cerr << "Error while releasing LP problem" << std::endl;
    }
  }
  status = CPXcloseCPLEX(&env);
  if( status ){
    std::cerr << "Error while releasing CPLEX environment" << std::endl;
    return 1;
  } 
  std::cout << "Library test ended successfully" << std::endl;
  return 0;
}

