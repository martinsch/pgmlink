/**
   @file
   @ingroup pgm
   @brief graphical model-based reasoner
*/

#ifndef REASONER_OPENGM_H
#define REASONER_OPENGM_H

#include <algorithm>
#include <iterator>
#include <map>
#include <vector>
#include <utility>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bimap.hpp>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include "pgmlink/event.h"
#include "pgmlink/feature.h"
#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"
#include "pgmlink/pgm_chaingraph.h"

namespace pgmlink {
  class Traxel;

  namespace pgm {
    typedef opengm::LPCplex<OpengmModel, opengm::Minimizer> OpengmLPCplex;
    using boost::shared_ptr;

    /**
       @brief Accessing entries of a Factor/Function that was already added to a graphical model.

       Manages a pointer to an element of an array-like opengm function (usually opengm::ExplicitFunction).
       Validity of the pointer is ensured by owning a smart pointer to the full model.

       Use this class to modify factor elements of an already instantiated opengm graphical model.
     */
    class FactorEntry {
    public:
    FactorEntry() : entry_(NULL) {}
    FactorEntry( shared_ptr<OpengmModel> m, /**< has to be valid */
		 OpengmModel::ValueType* entry /**< has to point into the opengm model to ensure the same lifetime */
		 ) :
      m_(m), entry_(entry) {}
      
      void set( OpengmModel::ValueType );
      OpengmModel::ValueType get() const;

      shared_ptr<OpengmModel> model() const { return m_; }

    private:
      shared_ptr<OpengmModel> m_;
      OpengmModel::ValueType* entry_;
    };
  } /* namespace pgm */


  class Chaingraph : public Reasoner {
    public:
    typedef pgm::chaingraph::Model::node_var_map node_var_map;
    typedef pgm::chaingraph::Model::arc_var_map arc_var_map;

    Chaingraph(bool with_constraints = true,
	       double ep_gap = 0.01,
	       bool fixed_detections = false
	       )
      : optimizer_(NULL),
      with_constraints_(with_constraints),
      fixed_detections_(fixed_detections),
      ep_gap_(ep_gap),
      builder_(NULL)
	{ builder_ = new pgm::chaingraph::ECCV12ModelBuilder(); (*builder_).with_detection_vars().with_divisions(); }
    

  Chaingraph(const pgm::chaingraph::ModelBuilder& builder,
	     bool with_constraints = true,
	     double ep_gap = 0.01,
	     bool fixed_detections = false
    ) 
    : optimizer_(NULL),
    with_constraints_(with_constraints),
    fixed_detections_(fixed_detections),
    ep_gap_(ep_gap),
    builder_(builder.clone())
    {};
    ~Chaingraph();

    virtual void formulate( const HypothesesGraph& );
    virtual void infer();
    virtual void conclude( HypothesesGraph& );

    double forbidden_cost() const;
    bool with_constraints() const;
    const pgm::chaingraph::ModelBuilder& builder() { return *builder_; }
    void builder(const pgm::chaingraph::ModelBuilder& builder) {
      if(builder_) delete builder_; builder_ = builder.clone(); }

    /** Return current state of graphical model
     *
     * The returned pointer may be NULL before formulate() is called
     * the first time.
     **/
    const pgm::OpengmModel* get_graphical_model() const;

    /** Return mapping from HypothesesGraph nodes to graphical model variable ids
     *
     * The map is populated after the first call to formulate().
     */
    const node_var_map& get_node_map() const;

    /** Return mapping from HypothesesGraph arcs to graphical model variable ids
     *
     * The map is populated after the first call to formulate().
     */
    const arc_var_map& get_arc_map() const;
    

    private:
    // copy and assingment have to be implemented, yet
    Chaingraph(const Chaingraph&) {};
    Chaingraph& operator=(const Chaingraph&) { return *this;};
    void reset();
    
    pgm::OpengmLPCplex* optimizer_;
    shared_ptr<pgm::chaingraph::Model> linking_model_;

    bool with_constraints_;
    bool fixed_detections_;

    double ep_gap_;
    pgm::chaingraph::ModelBuilder* builder_;
};
} /* namespace pgmlink */


/**/
/* implementation */
/**/
#include <boost/ptr_container/ptr_vector.hpp>

namespace pgmlink {
  template<class IT1, class IT2, class IT3>
    std::vector<pgm::OpengmModel::ValueType> pgm::chaingraph::ModelTrainer::train(IT1 samples_begin, IT1 samples_end, IT2 node_labels, IT3 arc_labels) const {
    // for each sample: build chaingraph model
    boost::ptr_vector<pgm::chaingraph::Model> models;
    chaingraph::ECCV12ModelBuilder b;
    b.with_detection_vars().with_divisions();
    for(IT1 sample=samples_begin; sample!=samples_end; ++sample){
      models.push_back(b.build(*sample));
    }

    // convert HypothesesGraph labels to OpengmModel labels
    IT2 cur_node_labels = node_labels;
    IT3 cur_arc_labels = arc_labels;
    std::vector<std::vector<pgm::OpengmModel::LabelType> > var_labels(std::distance(samples_begin, samples_end));
    for (int i=0; i < std::distance(samples_begin, samples_end); ++i) {
      var_labels[i] = std::vector<pgm::OpengmModel::LabelType>(models[i].opengm_model->numberOfVariables());
      for(pgm::chaingraph::Model::var_t var_idx = 0; var_idx < var_labels[i].size(); ++var_idx) {
	switch(models[i].var_category(var_idx)) {
	case pgm::chaingraph::Model::node_var: {
	  pgm::chaingraph::Model::node_t n = models[i].node_of_var(var_idx);
	  var_labels[i][var_idx] = (*cur_node_labels)[n];
	} break; 
	case pgm::chaingraph::Model::arc_var: {
	  pgm::chaingraph::Model::arc_t a = models[i].arc_of_var(var_idx);
	  var_labels[i][var_idx] = (*cur_arc_labels)[a];
	} break;
	  default:
	    throw std::runtime_error("chaingraph::ModelTrainer::train(): unknown var category encountered");
	    break;
	  }
      }
      ++cur_node_labels;
      ++cur_arc_labels;
    }
    
    return std::vector<pgm::OpengmModel::ValueType>();
  }

} /* namespace pgmlink */
#endif /* REASONER_OPENGM_H */
