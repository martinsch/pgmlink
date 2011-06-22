#include <cassert>
#include <exception>
#include <string>
#include "log.h"
#include "randomforest.h"

using namespace std;

namespace Tracking {
    namespace RF {

        vigra::RandomForest<RF_LABEL_TYPE> getRandomForest(std::string filename)
        {
            // check if file exists
            FILE* pFile = fopen ( filename.c_str(), "r" );
            if ( pFile == NULL){
                throw std::runtime_error("Tracking: Could not create Random Forest. Input file does not exist.");
            }

            // load the Random Forest
            vigra::RandomForest<RF_LABEL_TYPE> rf;
            vigra::rf_import_HDF5(rf, filename);

            return rf;
        }


        vigra::MultiArray<2,float> createFeatureVector(const Traxel &tr, const std::vector<std::string> &selFeatures)
        {
            // calculate the total size of the feature vector
            unsigned int len = 0;

            // go through all selected elements
            for(std::vector<std::string>::const_iterator it = selFeatures.begin(); it != selFeatures.end(); it++)
            {
                // find selected element in feature map
                FeatureMap::const_iterator f_it = tr.features.find(*it);
                if( f_it != tr.features.end() )
                {
                    // add its length, if found
                    len += (*f_it).second.size();
                }
                else
                {
                    // throw exception or print error message
                    //std::cerr << "Tracking: Feature '" << (*it) << "' not found for tracklet " << tr.Id << ".\n";
                }
            }

            // allocate memory for the feature matrix
            vigra::MultiArray<2,float> featureMatrix (matrix_shape(1,len));

            // fill matrix
            int featureCount = 0;
            for(std::vector<std::string>::const_iterator it = selFeatures.begin(); it != selFeatures.end(); it++)
            {
                // find selected element in feature map
                FeatureMap::const_iterator f_it = tr.features.find(*it);
                if( f_it != tr.features.end() )
                {
                    // copy entries
                    int len = (*f_it).second.size();
                    for(int i = 0; i < len; i++, featureCount++)
                    {
                        featureMatrix(0,featureCount) = (*f_it).second[i];
                    }
                }
                else
                {
                    std::cout << "Warning: Tracking::RF::createFeatureVector: feature "<< *it << " not found.\n";
                }
            }

            return featureMatrix;
        }
      
        vigra::MultiArray<2,double> getProbabilities(const vigra::MultiArray<2,float> &features, const vigra::RandomForest<RF_LABEL_TYPE> &rf)
        {
            // initialize probability matrix
            int len = features.shape(0);
            vigra::MultiArray<2,double> prob (matrix_shape(len,2));

            // predict probabilities
            rf.predictProbabilities(features, prob);

            return prob;
        }

      namespace {
	void save_as_feature(Traxel& tr, const string& name, feature_type value) {
	  Tracking::feature_array feat;
	  feat.push_back(value);
          tr.features[name] = feat;
	}
      }

      int predictTracklets( Traxels &ts,
                              vigra::RandomForest<RF_LABEL_TYPE> &rf,
                              std::vector<std::string> &sel,
                              unsigned int lbl = 1,
                              std::string lblname = "cellness") {
            int count = 0;
            for(Traxels::iterator tr = ts.begin(); tr != ts.end(); tr++, count++)
            {
	      double prob = predict(tr->second, rf, sel, lbl);
	      save_as_feature(tr->second, lblname, prob);
            }

            return count;
      }

      void predict_traxels( TraxelStore& ts,
                              const vigra::RandomForest<RF_LABEL_TYPE>& rf,
                              const std::vector<std::string>& feature_names,
                              unsigned int cls = 1,
                              const std::string& output_feat_name = "cellness") {
	for(TraxelStore::iterator it = ts.begin(); it != ts.end(); ++it) {
	  Traxel tr = *it;
	  double prob = predict(tr, rf, feature_names, cls);
	  save_as_feature(tr, output_feat_name, prob);
	  ts.replace(it, tr);
	}
      }

      double predict( const Traxel& tr, 
		      const vigra::RandomForest<RF_LABEL_TYPE>& rf, 
		      const std::vector<std::string>& feature_names,
		      unsigned int cls) {
	if(cls >= static_cast<unsigned int>(rf.class_count())) {
	  throw std::runtime_error("predict(): Provided class number is too large.");
	}

	// get the features of the current tracklet
	vigra::MultiArray<2,float> features ( createFeatureVector(tr, feature_names) );
	// predict the class probabilities
	vigra::MultiArray<2,double> prob = getProbabilities(features,rf);
	double ret = prob(0, cls);
	return ret;
      }




        /*
         * Helper function. Read a certain feature entry from file.
         */
        void loadFeature(vigra::HDF5File &f, Tracking::Traxels &ts, unsigned int label, std::string name, int size)
        {
            feature_array feature;
            feature.reshape(array_shape(size));
            std::stringstream path;
            path << "/features/" << label << "/" << name;

            try{
                f.read(path.str(),feature);
		Tracking::feature_array farr;
		for(feature_array::iterator it = feature.begin(); it != feature.end(); ++it) {
		  farr.push_back(*it);
		}
                ts[label].features[name] = farr;
                ts[label].Id = label;
            }
            catch(...){
                std::cerr << "Could not read feature '" << name << "' of label " << label <<".\n";
            }
        }

        Traxels loadTracklets(std::string filename)
        {
            // load features for every object & create feature matrix
            vigra::HDF5File f (filename, vigra::HDF5File::Open);
            int labelcount;
            f.read("/features/labelcount", labelcount);

            vigra::MultiArray<1,RF_LABEL_TYPE> labelcontent;
            labelcontent.reshape(array_shape(labelcount));
            f.read("/features/labelcontent",labelcontent);

            // init Traxels object
            Tracking::Traxels ts;

            // go through all objects and load features
            for(int i = 0; i < labelcount; i++){
                unsigned int label = i+1;
                if(labelcontent(i) == 1)
                {
                    // load normal object features (109)

                    //load feature volume
                    loadFeature(f,ts,label,"volume",1);

                    //load feature Bounding Box
                    loadFeature(f,ts,label,"bbox",7);

                    //load feature Unweighted Mean Position
                    loadFeature(f,ts,label,"position",12);

                    //load feature Weighted Mean Position
                    loadFeature(f,ts,label,"com",12);

                    //load feature Principal Component
                    loadFeature(f,ts,label,"pc",12);

                    //load feature Intensity
                    loadFeature(f,ts,label,"intensity",4);

                    //load feature Intensity Minimum/Maximum/Quantiles
                    loadFeature(f,ts,label,"intminmax",9);

                    //load feature Pairwise Energy
                    loadFeature(f,ts,label,"pair",4);

                    //load feature Statistical Geometric Features
                    loadFeature(f,ts,label,"sgf",48);

                    // load large features (89)

                    //load feature Weighted Mean Position
                    loadFeature(f,ts,label,"lcom",12);

                    //load feature Principal Component
                    loadFeature(f,ts,label,"lpc",12);

                    //load feature Intensity
                    loadFeature(f,ts,label,"lintensity",4);

                    //load feature Intensity Minimum/Maximum/Quantiles
                    loadFeature(f,ts,label,"lintminmax",9);

                    //load feature Pairwise Energy
                    loadFeature(f,ts,label,"lpair",4);

                    //load feature Statistical Geometric Features
                    loadFeature(f,ts,label,"lsgf",48);
                }
            }

            return ts;

        } /* loadTracklets */

    } /* namespace RF */
} /* namespace Tracking */
