/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# param_test.cpp
#
# Description: Test code for Param.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

namespace bdata = boost::unit_test::data;

#include "mpi.h"

#include "MolData.h"
#include "EmNNModel.h"
#include "LearnerTestMols.h"

#include <boost/filesystem.hpp>
#include <GraphMol/RDKitBase.h>

struct NNParamTestFixture {
    NNParamTestFixture() {

        std::vector<std::string> fnames;
        fnames.push_back("IonicFeatures");    //5 features + bias =  6 features
        // let us set a 6->4->2->1 NN
        // so there is 6 *4 + 5* 2 + 3 * 1= 24 + 10 + 3 = 37 weights
        std::vector<int> hlayer_numnodes = {4, 2, 1};
        std::vector<int> act_ids{RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION,
                                 RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION,
                                 LINEAR_NN_ACTIVATION_FUNCTION};

        std::vector<float> dropout_probs(2, 0.0);
        boost::container::vector<bool>is_frozen(3, false);
        param = new NNParam(fnames, 3, hlayer_numnodes, act_ids, dropout_probs, is_frozen);

        param->initWeights(PARAM_RANDOM_INIT);

        std::vector<float> weights(37, 0.0);

        // Input Layer 6 Nodes
        weights[0] = 1.0;
        weights[2] = -2.0;
        weights[3] = -4.0;
        weights[5] = 3.0;
        weights[6] = -3.0;
        weights[8] = 5.0;
        weights[9] = 12.0;
        weights[11] = -7.0;
        weights[12] = -2.0;
        weights[14] = 1.0;
        weights[15] = 6.0;
        weights[17] = -3.0;
        weights[18] = 5.0;
        weights[20] = 1.0;
        weights[21] = -2.0;
        weights[23] = 2.0;
        // Second Layer 4 Nodes + 1 bias Node
        weights[24] = 4.0;
        weights[25] = 2.0;
        weights[26] = 1.0;
        weights[27] = -3.0;
        weights[29] = -5.0;
        weights[30] = 1.0;
        weights[31] = 2.0;
        weights[32] = -2.0;
        // Second Layer 3 Nodes + 1 bias Node
        weights[34] = 5.0;
        weights[35] = -2.0;
        weights[36] = -1.0;    //Theta

        std::vector<float> weights_for_all_energy = weights;
        std::copy(weights.begin(), weights.end(), std::back_inserter(weights_for_all_energy));
        std::copy(weights.begin(), weights.end(), std::back_inserter(weights_for_all_energy));

        param->setWeights(weights_for_all_energy);

    };

    ~NNParamTestFixture() { delete param; };
    NNParam *param = nullptr;
};

BOOST_FIXTURE_TEST_SUITE(NNParamsComputeTests, NNParamTestFixture)

    BOOST_AUTO_TEST_CASE(ComputeThetasTest) {
        double tol = 1e-5;
        FeatureVector fv;

        fv.addFeatureAtIdx(1.0, 0);
        fv.addFeatureAtIdx(1.0, 2);
        fv.addFeatureAtIdx(1.0, 5);

        double theta = param->computeTheta(fv, 0);
        BOOST_CHECK_EQUAL(param->getNumWeightsPerEnergyLevel(), 37);
        BOOST_CHECK_EQUAL(param->getNumWeights(), 37 * 3);
        BOOST_CHECK_CLOSE_FRACTION(theta, 12.0, tol);
    }

    std::vector<int> energies = {0, 1, 2};

    BOOST_DATA_TEST_CASE(ComputeDeltasTest, bdata::make(energies), energy) {
        double tol = 1e-5;

        FeatureVector fv1, fv2;
        fv1.addFeatureAtIdx(1.0, 0);
        fv1.addFeatureAtIdx(1.0, 2);
        fv1.addFeatureAtIdx(1.0, 5);
        fv2.addFeatureAtIdx(1.0, 0);
        fv2.addFeatureAtIdx(1.0, 3);
        fv2.addFeatureAtIdx(1.0, 5);

        //Run the forwards paths
        std::vector<azd_vals_t> z_values(2), a_values(2);
        double theta1 = param->computeTheta(fv1, energy, z_values[0], a_values[0]);
        double theta2 = param->computeTheta(fv2, energy, z_values[1], a_values[1]);
        double rho_denom = 1.0 + exp(theta1) + exp(theta2);

        //Now run backwards to compute the delta values
        std::vector<azd_vals_t> deltasA, deltasB;
        param->computeDeltas(deltasA, deltasB, z_values, a_values, rho_denom, energy);

        //Check the values
        double x1 = exp(theta1) / rho_denom, x2 = exp(theta2) / rho_denom;
        double exp_deltaA_1[] = {-5.0, -4.0, 0.0, 0.0, -2.0, -1.0, 1.0};
        double exp_deltaB_1[] = {-5 * x1, -4 * x1, 0.0, 0.0, -2 * x1, -x1, x1};
        double exp_deltaA_2[] = {0.0, 0.0, 8.0, 0.0, -2.0, -1.0, 1.0};
        double exp_deltaB_2[] = {0.0, 0.0, 8 * x2, 0.0, -2 * x2, -x2, x2};

        BOOST_CHECK_EQUAL(deltasA.size(), 2);
        BOOST_CHECK_EQUAL(deltasB.size(), 2);
        BOOST_CHECK_EQUAL(deltasA[0].size(), 7);
        BOOST_CHECK_EQUAL(deltasB[0].size(), 7);
        BOOST_CHECK_EQUAL(deltasA[1].size(), 7);
        BOOST_CHECK_EQUAL(deltasB[1].size(), 7);

        for (int i = 0; i < 7; i++) {
            BOOST_CHECK_CLOSE_FRACTION(exp_deltaA_1[i], deltasA[0][i], tol);
            BOOST_CHECK_CLOSE_FRACTION(exp_deltaB_1[i], deltasB[0][i], tol);
            BOOST_CHECK_CLOSE_FRACTION(exp_deltaA_2[i], deltasA[1][i], tol);
            BOOST_CHECK_CLOSE_FRACTION(exp_deltaB_2[i], deltasB[1][i], tol);
        }
    }

    BOOST_DATA_TEST_CASE(UnweightedGradientsTest, bdata::make(energies), energy) {
        double tol = 1e-2;

        FeatureVector fv1, fv2;
        fv1.addFeatureAtIdx(1.0, 0);
        fv1.addFeatureAtIdx(1.0, 2);
        fv1.addFeatureAtIdx(1.0, 5);
        fv2.addFeatureAtIdx(1.0, 0);
        fv2.addFeatureAtIdx(1.0, 3);
        fv2.addFeatureAtIdx(1.0, 5);

        //Run the forwards paths
        std::vector<azd_vals_t> z_values(2), a_values(2);
        double theta1 = param->computeTheta(fv1, energy, z_values[0], a_values[0]);
        double theta2 = param->computeTheta(fv2, energy, z_values[1], a_values[1]);
        double rho_denom = 1.0 + exp(theta1) + exp(theta2);

        //Now run backwards to compute the delta values
        std::vector<azd_vals_t> deltasA, deltasB;
        param->computeDeltas(deltasA, deltasB, z_values, a_values, rho_denom, energy);

        //Now compute the unweighted gradients
        std::vector<std::vector<float> > unweighted_grads;
        std::set<unsigned int> used_idxs;
        std::vector<const FeatureVector *> fvs(2);
        fvs[0] = &fv1;
        fvs[1] = &fv2;
        param->computeUnweightedGradients(unweighted_grads, used_idxs, fvs, deltasA, deltasB, a_values);

        //Check the values
        float x1 = exp(theta1) / rho_denom, x2 = exp(theta2) / rho_denom;
        float exp_unweighted_1[] = {-5 + 5 * x1, 0, -5 + 5 * x1, 0, 0, -5 + 5 * x1, -4 + 4 * x1, 0, -4 + 4 * x1, 0, 0,
                                     -4 + 4 * x1, -8 * x2, 0, 0, -8 * x2, 0, -8 * x2, 0, 0, 0, 0, 0, 0,
                                     -2 + 2 * x1 + 2 * x2, -2 * 2 + 2 * 2 * x1, 5 * 2 - 5 * 2 * x1, 2 * x2, 0,
                                     -1 + x1 + x2, -2 + 2 * x1, 5 - 5 * x1, x2, 0, 1 - x1 - x2, 3 - 3 * x1 - x2,
                                     -13 + 13 * x1 + 7 * x2};
        float exp_unweighted_2[] = {5 * x1, 0, 5 * x1, 0, 0, 5 * x1, 4 * x1, 0, 4 * x1, 0, 0, 4 * x1, 8 - 8 * x2, 0, 0,
                                     8 - 8 * x2, 0, 8 - 8 * x2, 0, 0, 0, 0, 0, 0, -2 + 2 * x1 + 2 * x2, 4 * x1,
                                     -5 * 2 * x1, -2 + 2 * x2, 0, -1 + x1 + x2, 2 * x1, -5 * x1, -1 + x2, 0,
                                     1 - x1 - x2, +1 - 3 * x1 - x2, -7 + 13 * x1 + 7 * x2};
        float exp_unweighted_persist[] = {5 * x1, 0, 5 * x1, 0, 0, 5 * x1, 4 * x1, 0, 4 * x1, 0, 0, 4 * x1, -8 * x2, 0,
                                           0, -8 * x2, 0, -8 * x2, 0, 0, 0, 0, 0, 0, 2 * x1 + 2 * x2, 4 * x1,
                                           -5 * 2 * x1, 2 * x2, 0, x1 + x2, 2 * x1, -5 * x1, x2, 0, -x1 - x2,
                                           -3 * x1 - x2, 13 * x1 + 7 * x2};

        BOOST_CHECK_EQUAL(unweighted_grads.size(), 3);
        BOOST_CHECK_EQUAL(unweighted_grads[0].size(), 37);
        BOOST_CHECK_EQUAL(unweighted_grads[1].size(), 37);
        BOOST_CHECK_EQUAL(unweighted_grads[2].size(), 37);
        BOOST_CHECK_EQUAL(used_idxs.size(), 16);

        for (int i = 0; i < 37; i++) {
            BOOST_CHECK_CLOSE_FRACTION(exp_unweighted_1[i], unweighted_grads[0][i], tol);
            BOOST_CHECK_CLOSE_FRACTION(exp_unweighted_2[i], unweighted_grads[1][i], tol);
            BOOST_CHECK_CLOSE_FRACTION(exp_unweighted_persist[i], unweighted_grads[2][i], tol);
        }

        std::set<int> expected_used_idxs = {0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23};
        BOOST_TEST(expected_used_idxs == used_idxs, boost::test_tools::per_element());
    }

    BOOST_AUTO_TEST_CASE(GradientTest) {
        // set up mpi
        int mp_init_flag;
        MPI_Initialized(&mp_init_flag);
        if(mp_init_flag == 0)
            MPI::Init();

        double tol = 1e-3;
        std::vector<float> grads;

        suft_counts_t suft;

        //Model Config
        config_t cfg; initDefaultConfig(cfg);
        cfg.model_depth = 2;
        cfg.spectrum_depths.push_back(2);
        cfg.spectrum_weights.push_back(1);
        cfg.lambda = 0.01;
        initDerivedConfig(cfg);

        //Create the molecule data
        NNParamTestMol moldata(&cfg);

        //Set some parameter weights
        std::vector<std::string> fnames;
        fnames.push_back("IonicFeatures");

        std::string param_filename = "tmp_param_file.log";
        if( boost::filesystem::exists( param_filename ) )
            boost::filesystem::remove( param_filename );
        param->saveToFile( param_filename );

        //Set some arbitrary suft values
        suft.values.resize(1);

        unsigned int N = moldata.getNumTransitions() + moldata.getNumFragments();
        suft.values[0].assign(3*N, 0.0);
        suft.values[0][0] = 0.2; suft.values[0][1] = 0.3; suft.values[0][2] = 0.5; suft.values[0][3] = 0.0; suft.values[0][4] = 0.0;
        //suft.values[0][N] = 0.9; suft.values[0][N+1] = 0.05; suft.values[0][N+2] = 0.05; suft.values[0][N+3] = 0.0; suft.values[0][N+4] = 0.0;
        //suft.values[0][2*N] = 0.08; suft.values[0][2*N+1] = 0.9; suft.values[0][2*N+2] = 0.02; suft.values[0][2*N+3] = 0.0; suft.values[0][2*N+4] = 0.0;

        //Initialise all gradients to 0
        grads.resize(param->getNumWeights());
        for( unsigned int i = 0; i < grads.size(); i++ )
            grads[i] = 0.0;

        std::set<unsigned int> used_idxs;
        FeatureCalculator fc_null(fnames);
        std::string null_str = "null";
        EmNNModel em(&cfg, &fc_null, null_str, param_filename );
        int energy_level = 0;
        int mol_idx = 0;
        double Q_only  = em.computeLogLikelihoodLoss(mol_idx, moldata, suft, energy_level);
        em.collectUsedIdx(moldata,used_idxs, 0);
        em.computeAndAccumulateGradient(&grads[0], 0, moldata, suft, 0, energy_level);

        //Check Q
        double theta1 = 12.0, theta2 = 10.0;
        double rho_denom = 1.0 + exp(theta1) + exp(theta2);
        double x1 = exp(theta1)/rho_denom, x2 = exp(theta2)/rho_denom;

        double expected_Q = (0.2)*(theta1-log(rho_denom)) + (0.3)*(theta2-log(rho_denom)) + (0.5)*(-log(rho_denom));

        //Check the gradients
        double exp_unweighted_1[] = {-5+5*x1,0,-5+5*x1,0,0,-5+5*x1,-4+4*x1,0,-4+4*x1,0,0,-4+4*x1,-8*x2, 0, 0, -8*x2,0,-8*x2,0,0,0,0,0,0,-2+2*x1+2*x2,-2*2+2*2*x1, 5*2-5*2*x1,2*x2,0,-1+x1+x2,-2+2*x1,5-5*x1,x2,0,1-x1-x2,3-3*x1-x2, -13+13*x1+7*x2 };
        double exp_unweighted_2[] = {5*x1,0,5*x1,0,0,5*x1, 4*x1,0,4*x1,0,0,4*x1, 8-8*x2, 0,0,8-8*x2,0,8-8*x2, 0,0,0,0,0,0, -2+2*x1+2*x2,4*x1,-5*2*x1,-2+2*x2,0,-1+x1+x2, 2*x1,-5*x1, -1+x2, 0, 1-x1-x2,+1-3*x1-x2,-7+13*x1+7*x2 };
        double exp_unweighted_persist[] = {5*x1,0,5*x1,0,0,5*x1, 4*x1,0,4*x1,0,0,4*x1,-8*x2, 0,0,-8*x2,0,-8*x2,0,0,0,0,0,0,2*x1+2*x2,4*x1,-5*2*x1,2*x2,0,x1+x2, 2*x1,-5*x1, x2, 0, -x1-x2,-3*x1-x2,13*x1+7*x2 };


        BOOST_CHECK_CLOSE_FRACTION(expected_Q, Q_only, tol);

        //Energy 0
        for( unsigned int i = 0; i < param->getNumWeightsPerEnergyLevel(); i++ )
            BOOST_CHECK_CLOSE_FRACTION(grads[i], (0.2*exp_unweighted_1[i] + 0.3*exp_unweighted_2[i] + 0.5*exp_unweighted_persist[i]), tol);

        //Check used Index
        std::vector<unsigned int> expected_idxs = {0,2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36};
        BOOST_CHECK_EQUAL_COLLECTIONS(expected_idxs.begin(), expected_idxs.end(),
                                      used_idxs.begin(), used_idxs.end());

    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(NNParamsIndexAndDropOutTests)

    BOOST_AUTO_TEST_CASE(NNParamsTestBiasIndexes) {
        std::vector<std::string> fnames;
        fnames.push_back("IonicFeatures");    //5 features + bias =  6 features
        std::vector<int> hlayer_numnodes = {4, 2, 2, 1};
        std::vector<int> act_ids{RELU_NN_ACTIVATION_FUNCTION, RELU_NN_ACTIVATION_FUNCTION,
                                 RELU_NN_ACTIVATION_FUNCTION, LINEAR_NN_ACTIVATION_FUNCTION};

        std::vector<float> dropout_probs(2, 0.0);
        NNParam param(fnames, 1, hlayer_numnodes, act_ids, dropout_probs);

        std::vector<unsigned int> bias_indexes;
        param.getBiasIndexes(bias_indexes);
        BOOST_CHECK_EQUAL(bias_indexes.size(), 9);

        std::vector<unsigned int> expected_idxs = {0, 6, 12, 18, 24, 29, 34, 37, 40};
        BOOST_CHECK_EQUAL_COLLECTIONS(expected_idxs.begin(), expected_idxs.end(),
                                      bias_indexes.begin(), bias_indexes.end());
    }

    BOOST_AUTO_TEST_CASE(NNParamsTestDropout) {
        std::vector<std::string> fnames;
        fnames.push_back("IonicFeatures");    //5 features + bias =  6 features
        std::vector<int> hlayer_numnodes{64, 32, 32, 1};
        std::vector<int> act_ids{RELU_NN_ACTIVATION_FUNCTION, RELU_NN_ACTIVATION_FUNCTION,
                                 RELU_NN_ACTIVATION_FUNCTION, LINEAR_NN_ACTIVATION_FUNCTION};

        std::vector<float> dropout_probs{0.5, 0.5, 0, 0};
        NNParam param(fnames, 1, hlayer_numnodes, act_ids, dropout_probs);

        auto dropout_prob_ptr = param.getDropoutsProbPtr();
        BOOST_CHECK_EQUAL(dropout_prob_ptr->size(), dropout_probs.size());

        BOOST_CHECK_EQUAL_COLLECTIONS(dropout_probs.begin(), dropout_probs.end(),
                                      dropout_prob_ptr->begin(), dropout_prob_ptr->end());

        std::string prev_selection = "";
        std::string current_selection = "";

        for (int trail = 0; trail < 10; ++trail) {
            param.rollDropouts();
            auto dropout_ptr = param.getDropoutsPtr();
            auto neuron_idx = 0;
            prev_selection = current_selection;

            for (auto h_layer_idx = 0; h_layer_idx < hlayer_numnodes.size(); ++h_layer_idx) {
                int num_active_neuron = 0;
                auto num_neuron = hlayer_numnodes[h_layer_idx];
                for (int i = 0; i < num_neuron; ++i) {
                    if (!(*dropout_ptr)[neuron_idx]) {
                        num_active_neuron++;
                        current_selection += '0';
                    } else {
                        current_selection += '1';
                    }
                    neuron_idx++;
                }

                double expected_active_ratio = 1.0 - dropout_probs[h_layer_idx];
                double active_ratio = (double) num_active_neuron / (double) num_neuron;
                BOOST_CHECK_CLOSE_FRACTION(expected_active_ratio, active_ratio, 0.05);
                BOOST_CHECK(current_selection != prev_selection);
            }
        }
    }
    BOOST_AUTO_TEST_CASE(SaveAndLoadParamFileTest){

        double tol = 1e-5;

        //Create some parameters
        std::vector<std::string> fnames;
        fnames.push_back("IonicFeatures"); // 6 features in total

        std::vector<int> hlayer_numnodes = {12,4,1};
        std::vector<int> act_ids{RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION,
                                 RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION,
                                 LINEAR_NN_ACTIVATION_FUNCTION};

        std::vector<float> dropout_probs(2,0.0);
        NNParam param(fnames, 1, hlayer_numnodes, act_ids, dropout_probs);
        param.initWeights(PARAM_RANDOM_INIT);

        //Save to file
        std::string filename = "tmp_param_file.log";
        if( boost::filesystem::exists( filename ) )
            boost::filesystem::remove( filename );
        param.saveToFile( filename );

        //Load from file
        NNParam param_load( filename );

        //Check the configuration
        BOOST_CHECK_EQUAL(param.getNumEnergyLevels(), param_load.getNumEnergyLevels());

        //Check the feature names
        std::vector<std::string> *fnames1 = param.getFeatureNames();
        std::vector<std::string> *fnames2 = param_load.getFeatureNames();
        BOOST_TEST((*fnames1) == (*fnames2 ), boost::test_tools::per_element());

        //Check the weights
        BOOST_CHECK_EQUAL(param.getNumWeights(), param_load.getNumWeights());

        for( int i = 0; i < param.getNumWeights(); i++ )
            BOOST_CHECK_CLOSE_FRACTION(param.getWeightAtIdx(i),param_load.getWeightAtIdx(i),tol);

        //Check that a computed theta value is the same
        FeatureVector fv;
        fv.addFeatureAtIdx(1.0, 0); fv.addFeatureAtIdx(1.0, 2); fv.addFeatureAtIdx(1.0, 5);

        double theta1 = param.computeTheta( fv, 0 );
        double theta2 = param_load.computeTheta( fv, 0 );
        BOOST_CHECK_CLOSE_FRACTION(theta1,theta2,tol);

    }
BOOST_AUTO_TEST_SUITE_END()