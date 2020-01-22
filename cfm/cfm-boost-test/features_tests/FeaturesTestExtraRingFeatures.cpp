/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for features.cpp
#
# Author: Felicity Allen, Fei Wang
# Created: November 2012, 2019
#########################################################################*/
#include "FeaturesTestsUtil.h"

BOOST_AUTO_TEST_SUITE(FeaturesTestExtraRingFeatures)

    BOOST_AUTO_TEST_CASE(BetweenTwoRings) {

        std::string smiles_or_inchi("c1ccccc1-C2CCCCC2-CCC");
        int break_id = 0;

        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "ExtraRingFeatures", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 1, 2, 3};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(SingleBondBetweenRingAndNonRing) {
        std::string smiles_or_inchi("c1ccccc1-C2CCCCC2-CCC");
        int break_id = 1;

        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "ExtraRingFeatures", break_id, false);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 1, 3};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(BreakRing) {
        //No Test Case Yet
        //TODO: Fix Me
    }

BOOST_AUTO_TEST_SUITE_END()