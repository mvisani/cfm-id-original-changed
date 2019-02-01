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
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "mpi.h"

#include "MolData.h"
#include "EmModel.h"
#include "LearnerTestMols.h"

#include <boost/filesystem.hpp>
#include <GraphMol/RDKitBase.h>

struct LinearModelThetaTestFixture {
    LinearModelThetaTestFixture() {

        std::vector<std::string> fnames;
        fnames.push_back("BreakAtomPair");
        param = new Param(fnames, 3);

        int med_offset = param->getNumWeightsPerEnergyLevel();

        //Set some parameter weights
        param->setWeightAtIdx(0.234, 0);
        param->setWeightAtIdx(-3.124, 1);
        param->setWeightAtIdx(0.89, 2);

        param->setWeightAtIdx(-5.23432, 0 + med_offset);
        param->setWeightAtIdx(0.0, 1 + med_offset);
        param->setWeightAtIdx(1.0, 2 + med_offset);

        int high_offset = 2 * med_offset;
        param->setWeightAtIdx(-35.0, 0 + high_offset);
        param->setWeightAtIdx(-10.0, 1 + high_offset);
        param->setWeightAtIdx(0.0, 2 + high_offset);
    };

    ~LinearModelThetaTestFixture() { delete param; };
    Param *param = nullptr;
};

BOOST_FIXTURE_TEST_SUITE(ParamsTestComputeTransitionThetas, LinearModelThetaTestFixture)

    BOOST_AUTO_TEST_CASE(FeatureVectorOne) {
        double tol = 1e-10;
        FeatureVector fv;

        int med_offset = param->getNumWeightsPerEnergyLevel();

        fv.addFeatureAtIdx(0.0, med_offset - 1);    //Set the feature length
        fv.addFeatureAtIdx(1.0, 1);
        fv.addFeatureAtIdx(1.0, 2);

        //Compute the transition thetas
        double theta_low = param->computeTheta(fv, 0);
        double theta_med = param->computeTheta(fv, 1);
        double theta_high = param->computeTheta(fv, 2);

        BOOST_CHECK_CLOSE_FRACTION(theta_low, -2.234, tol);
        BOOST_CHECK_CLOSE_FRACTION(theta_med, 1.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(theta_high, -10.0, tol);

    }

    BOOST_AUTO_TEST_CASE(FeatureVectorTwo) {
        double tol = 1e-10;
        FeatureVector fv;

        int med_offset = param->getNumWeightsPerEnergyLevel();

        fv.addFeatureAtIdx(0.0, med_offset - 1);    //Set the feature length
        fv.addFeatureAtIdx(1.0, 0);


        //Compute the transition thetas
        double theta_low = param->computeTheta(fv, 0);
        double theta_med = param->computeTheta(fv, 1);
        double theta_high = param->computeTheta(fv, 2);

        BOOST_CHECK_CLOSE_FRACTION(theta_low, 0.234, tol);
        BOOST_CHECK_CLOSE_FRACTION(theta_med, -5.23432, tol);
        BOOST_CHECK_CLOSE_FRACTION(theta_high, -35.0, tol);

    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ParamsTestComputeAndAccumulateGradient)

    BOOST_AUTO_TEST_CASE(GradientTest) {
        // set up mpi
        int mp_init_flag;
        MPI_Initialized(&mp_init_flag);
        if(mp_init_flag == 0)
            MPI::Init();

        double tol = 1e-3;

        std::vector<double> grads;

        suft_counts_t suft;

        //Model Config
        config_t cfg;
        cfg.model_depth = 3;
        int tmp_array[3] = {1, 2, 3};
        cfg.spectrum_depths.assign(tmp_array, tmp_array + 3);
        cfg.lambda = 0.01;
        initDerivedConfig(cfg);

        //Create the molecule data
        ParamTestMol moldata(&cfg);
        //Set some arbitrary parameter weights (only used in the regularizer term and for sizing)
        std::vector<std::string> fnames;
        fnames.push_back("HydrogenMovement");    //Size = 10
        Param param(fnames, 1);
        for (unsigned int i = 0; i < param.getNumWeights(); i++)
            param.setWeightAtIdx(0.0, i);

        param.setWeightAtIdx(0.2, 0);
        param.setWeightAtIdx(-0.5, 1);

        std::string param_filename = "tmp_param_file.log";
        if (boost::filesystem::exists(param_filename))
            boost::filesystem::remove(param_filename);
        param.saveToFile(param_filename);

        //Set some arbitrary suft values
        suft.values.resize(1);

        unsigned int N = moldata.getNumTransitions() + moldata.getNumFragments();
        suft.values[0].assign(N, 0.0);
        suft.values[0][1] = 0.5;
        suft.values[0][7] = 0.5;

        //Initialise all gradients to 0
        grads.resize(11);
        std::fill(grads.begin(), grads.end(), 0);


        std::set<unsigned int> used_idxs;
        FeatureCalculator fc_null(fnames);
        std::string null_str = "null";
        EmModel em(&cfg, &fc_null, null_str, param_filename);

        //Check Q
        double Q = em.computeLogLikelihoodLoss(0, moldata, suft, 0);
        double expected_Q = -1.236;
        BOOST_CHECK_CLOSE_FRACTION(Q, expected_Q, tol);

        // get used flags
        em.computeAndAccumulateGradient(&grads[0], 0, moldata, suft, true, used_idxs, 0, 0);
        //Check the used flags
        BOOST_CHECK_EQUAL(used_idxs.find(0) != used_idxs.end(), true);
        BOOST_CHECK_EQUAL(used_idxs.find(1) != used_idxs.end(), true);
        //Check one that shouldn't be on
        BOOST_CHECK_EQUAL(used_idxs.find(5) == used_idxs.end(), true);

        //Check the gradients for low energy level
        int energy = 0;
        em.computeAndAccumulateGradient(&grads[0], 0, moldata, suft, false, used_idxs, 0, energy);
        double expected_grads[6] = {-0.1624, 0.2499, 0, 0, 0, 0}; //-0.0568,0.2554,0.1649,0.4663};

        for (unsigned int i = 0; i < 6; i++) {
            BOOST_CHECK_CLOSE_FRACTION(grads[i], expected_grads[i], tol);
        }

    }

BOOST_AUTO_TEST_SUITE_END()
