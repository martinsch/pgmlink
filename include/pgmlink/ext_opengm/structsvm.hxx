#pragma once
#ifndef OPENGM_STRUCTSVM_HXX
#define OPENGM_STRUCTSVM_HXX

#include <cassert>
#include <ostream>
#include <vector>

#include <dlib/svm.h>
#include <dlib/optimization.h>
#include <dlib/matrix.h>

#include "opengm/opengm.hxx"
#include "pgmlink/ext_opengm/loss_hamming.hxx"
#include "pgmlink/ext_opengm/loglinearmodel.hxx"

#ifdef WITH_GUROBI
#include "opengm/inference/lpgurobi.hxx"
#else
#include "opengm/inference/lpcplex.hxx"
#endif

#include "opengm/operations/minimizer.hxx"


namespace opengm {

template<class LLM>
class StructSvmDlib
  : public dlib::structural_svm_problem<dlib::matrix<typename LLM::ValueType> > {
public:
  typedef LLM LoglinearModelType;

  StructSvmDlib(const std::vector<LoglinearModelType>& samples,
		const std::vector<std::vector<typename LoglinearModelType::LabelType> > labels)
    : samples_(samples), samples_and_loss_(samples), labels_(labels) {
    for(size_t i=0; i < samples.size(); ++i) {
      augmentWithLoss(samples_and_loss_[i], labels_[i]);
      LLM loss(samples_[i].space());
      augmentWithLoss(loss, labels_[i]);
      loss_.push_back(loss);
    }
  };

  void train(std::vector<typename LoglinearModelType::ValueType>& learned_weights);

  ///
  /// interface required by dlib
  ///
  typedef dlib::matrix<typename LLM::ValueType> matrix_type_dlib;
  typedef typename matrix_type_dlib::type scalar_type_dlib;
  typedef matrix_type_dlib feature_vector_type_dlib;

  virtual long get_num_dimensions () const;
  virtual long get_num_samples () const ;
  virtual void get_truth_joint_feature_vector ( long idx, feature_vector_type_dlib& psi ) const;
  virtual void separation_oracle (
				  const long idx,
				  const matrix_type_dlib& current_solution,
				  scalar_type_dlib& loss,
				  feature_vector_type_dlib& psi
				  ) const;

private:
  void augmentWithLoss( LoglinearModelType&, std::vector<typename LoglinearModelType::LabelType>& );

  mutable std::vector<LoglinearModelType> samples_;
  mutable std::vector<LoglinearModelType> samples_and_loss_;
  std::vector<LoglinearModelType> loss_;
  std::vector<std::vector<typename LoglinearModelType::LabelType> > labels_;
};
  
/******************/
/* implementation */
/******************/
  
template<class LLM>
void StructSvmDlib<LLM>::train(std::vector<typename LoglinearModelType::ValueType>& learned_weights) {
  learned_weights.clear();
  dlib::matrix<typename LLM::ValueType> weights;
  dlib::oca solver;
  solver.set_subproblem_epsilon(0.0000001);
  //solver.set_subproblem_max_iterations(100);
  std::cout << "\n";
  std::cout << "oca: subproblem epsilon: " << solver.get_subproblem_epsilon() << "\n";
  std::cout << "oca: subproblem max iterations: " << solver.get_subproblem_max_iterations() << "\n";

  this->be_verbose();
  solver( *this, weights);
  //assert(weights.size() == 2);
  assert(weights.nc() == 1);
  for(long i = 0; i< weights.nr(); ++i) {
    learned_weights.push_back(weights(0,i));
  }
}

template<class LLM>
long StructSvmDlib<LLM>::get_num_dimensions() const {
  return samples_.front().numberOfWeights();
}

template<class LLM>
long StructSvmDlib<LLM>::get_num_samples() const {
  return samples_.size();
}

template<class LLM>
void StructSvmDlib<LLM>::get_truth_joint_feature_vector( long idx, feature_vector_type_dlib& psi ) const {
  psi.set_size(this->get_num_dimensions(), 1);

  std::vector<typename LLM::ValueType> feats(samples_[idx].numberOfWeights());
  samples_[idx].weightedFeatureSums( labels_[idx], feats);
  for(long i = 0; i < psi.size(); ++i) {
    psi(i,0) = feats[i];
    std::cout << "true feature["<<i<<"]: " << psi(i,0) << "\n";
  }
}

template<class LLM>
void StructSvmDlib<LLM>::separation_oracle(const long idx,
					   const matrix_type_dlib& current_solution,
					   scalar_type_dlib& loss,
					   feature_vector_type_dlib& psi
					   ) const {
  std::cout << "\n";
  std::cout << "Entered separation oracle\n";
  std::vector<typename LLM::ValueType> weights;
  for(int i = 0; i < current_solution.nr(); ++i) {
    weights.push_back(-1*current_solution(i, 0));
    std::cout << "weights["<<i<<"]: " << current_solution(i,0) << "\n";
  }
  samples_and_loss_[idx].setWeights(weights);
  samples_[idx].setWeights(weights);

#ifdef WITH_GUROBI
  typedef opengm::LPGurobi<LLM, opengm::Minimizer> optimizer_type;
#else
  typedef opengm::LPCplex<LLM, opengm::Minimizer> optimizer_type;
#endif
  typename optimizer_type::Parameter params;
  params.verbose_ = true;
  params.integerConstraint_ = true;
  params.epGap_ = 0.01;
  opengm::InferenceTermination status;  
  optimizer_type inference = optimizer_type(samples_and_loss_[idx], params);
  status = inference.infer();
  if(status != opengm::NORMAL) {
    throw RuntimeError("GraphicalModel::infer(): optimizer terminated unnormally");
  }


  std::vector<typename LLM::LabelType> optimal;
  status = inference.arg(optimal);
  if(status != opengm::NORMAL) {
    throw RuntimeError("GraphicalModel::infer(): solution extraction terminated unnormally");
  }
  std::cout << "current[0]: " << optimal[0] << " label[0]: " << labels_[idx][0]  << "\n";
  std::cout << "current[1]: " << optimal[1] << " label[1]: " << labels_[idx][1]  << "\n";
  std::cout << "current[2]: " << optimal[2] << " label[2]: " << labels_[idx][2]  << "\n";
  std::cout << "current[3]: " << optimal[3] << " label[3]: " << labels_[idx][3]  << "\n";
  std::cout << "etc. etc.\n"; 

  std::vector<typename LLM::ValueType> feats(samples_and_loss_[idx].numberOfWeights());
  samples_and_loss_[idx].weightedFeatureSums( optimal, feats );
  psi.set_size(samples_[idx].numberOfWeights(), 1);
  for(long i = 0; i < psi.size(); ++i) {
    std::cout << "feature["<<i <<"]: " << feats[i] << "\n";
    psi(i,0) = feats[i];
  }

  loss = -1*loss_[idx].evaluate(optimal);
  std::cout << "ground truth energy for idx" << idx <<": "<< -1*samples_[idx].evaluate(labels_[idx])<<"\n";
  std::cout << "energy for idx" << idx <<": "<< -1*samples_[idx].evaluate(optimal)<<"\n";
  std::cout << "energy incl. loss for idx" << idx <<": "<< -1*samples_and_loss_[idx].evaluate(optimal)<<"\n";
  std::cout << "loss for idx "<< idx <<": " << loss << "\n";
  std::cout << "leaving oracle\n\n";
}

template<class LLM>
void StructSvmDlib<LLM>::augmentWithLoss( LLM& m, std::vector<typename LoglinearModelType::LabelType>& target_labels ) {
  typename LLM::IndexType lidx[] = {0};
  HammingFunction<typename LLM::ValueType, typename LLM::IndexType, typename LLM::LabelType> hf;
  for(typename LLM::IndexType i = 0; i < target_labels.size(); ++i) {
    lidx[0] = i;
    hf = HammingFunction<typename LLM::ValueType, typename LLM::IndexType, typename LLM::LabelType>(m.numberOfLabels(i), target_labels[i]);
    m.addFactor(m.addFunction(hf), lidx, lidx+1 );
  }
}

} // end namespace opengm

#endif
