/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for param->cpp
#
# Author: Felicity Allen, Fei Wang
# Created: November 2012, 2019
#########################################################################*/
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Param.h"
#include "Solver.h"

#include <GraphMol/RDKitBase.h>
#include <boost/filesystem.hpp>

struct SolverTestFixture {
  SolverTestFixture() {
    grads = {1.0, 0.0, -3.0, 3.3, -1.0, 0.0, 20, 17.8};
    used_idxs = {0, 1, 2, 3};

    std::vector<float> init_weights(8, 0.0);
    std::vector<std::string> fake_fnames = {"BreakAtomPair"};
    int num_energy_levels = 1;
    param = boost::shared_ptr<Param>(new Param(fake_fnames, num_energy_levels));
    param->setWeights(init_weights);
  };
  ~SolverTestFixture(){};
  std::vector<float> grads;
  boost::shared_ptr<Param> param;
  std::set<unsigned int> used_idxs;
};

BOOST_FIXTURE_TEST_SUITE(SolverTests, SolverTestFixture)

BOOST_AUTO_TEST_CASE(SGDTests) {
  double learning_rate = 0.1;
  double tol = 1e-6;
  Sgd Solver(learning_rate);
  // only use first half
  Solver.adjustWeights(grads, used_idxs, param);
  auto weight_ptr = param->getWeightsPtr()->begin();
  std::vector<double> expected_weights = {0.1, 0.0, -0.3, 0.33,
                                          0.0, 0.0, 0.0,  0.0};
  for (int i = 0; i < expected_weights.size(); ++i) {
    double w = *weight_ptr;
    BOOST_CHECK_CLOSE_FRACTION(w, expected_weights[i], tol);
    weight_ptr++;
  }
}

BOOST_AUTO_TEST_CASE(MomentumTests) {
  double tol = 1e-6;

  double learning_rate = 0.1;
  double momentum = 0.9;
  Momentum Solver(grads.size(), learning_rate, momentum);

  // iteration one, when there momentum doing nothing
  Solver.adjustWeights(grads, used_idxs, param);
  auto weight_ptr = param->getWeightsPtr()->begin();
  std::vector<double> expected_weights = {0.1, 0.0, -0.3, 0.33,
                                          0.0, 0.0, 0.0,  0.0};
  for (int i = 0; i < expected_weights.size(); ++i) {
    double w = *weight_ptr;
    BOOST_CHECK_CLOSE_FRACTION(w, expected_weights[i], tol);
    weight_ptr++;
  }

  // iteration one, when there momentum doing nothing
  Solver.adjustWeights(grads, used_idxs, param);
  weight_ptr = param->getWeightsPtr()->begin();
  double factor = 1 + 1 + momentum;
  expected_weights = {0.1 * factor, 0.0, -0.3 * factor, 0.33 * factor,
                      0.0,          0.0, 0.0,           0.0};
  for (int i = 0; i < expected_weights.size(); ++i) {
    double w = *weight_ptr;
    BOOST_CHECK_CLOSE_FRACTION(w, expected_weights[i], tol);
    weight_ptr++;
  }
}

BOOST_AUTO_TEST_CASE(AdamTests) {
  double tol = 1e-6;

  double learning_rate = 0.1;
  double beta_1 = 0.9;
  double beta_2 = 0.9;
  double eps = 1e-6;
  // Adam Solver (grads.size(), learning_rate, beta_1, beta_2, eps);

  // iteration one, when there momentum doing nothing
  // m_t = 0.1 * grads
  // v_t = 0.1 * grads^2
  // hat_m = m_t / ( 1 - beta_1 ^ t) = grads
  // hat_v= v_t / ( 1 - beta_2 ^ t) = grads
  // w = lr * m_hat / ( sqrt(v_hat) + eps) = grads / (sqrt(grad) + eps) * lr
  /* Solver.adjustWeights(grads, used_idxs, param);
   auto weight_ptr = param->getWeightsPtr()->begin();
   std::vector<double> expected_weights = { 0.1, 0.0, -0.3, 0.33, 0.0, 0.0, 0.0,
   0.0}; for(int i = 0; i < expected_weights.size(); ++i){ double w =
   *weight_ptr; BOOST_CHECK_CLOSE_FRACTION(w, expected_weights[i], tol);
       weight_ptr++;
   }*/

  // iteration one, when there momentum doing nothing
  /*Solver.adjustWeights(grads, used_idxs, param);
  weight_ptr = param->getWeightsPtr()->begin();
  double factor = 1 + 1 + momentum;
  expected_weights = { 0.1 * factor , 0.0, -0.3 * factor, 0.33 * factor, 0.0,
  0.0, 0.0, 0.0}; for(int i = 0; i < expected_weights.size(); ++i){ double w =
  *weight_ptr; BOOST_CHECK_CLOSE_FRACTION(w, expected_weights[i], tol);
      weight_ptr++;
  }*/
}

BOOST_AUTO_TEST_SUITE_END()