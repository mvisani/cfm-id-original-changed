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

BOOST_AUTO_TEST_SUITE(FeaturesTestRadicalFeatures)

    BOOST_AUTO_TEST_CASE(IonRadicalOnly) {
        FeatureVector *fv = getFeatureVector("CC[CH3+]", "CC(=O)O", "RadicalFeatures");
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 1};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(NLRadicalOnly) {
        FeatureVector *fv = getFeatureVector("CCC[CH4+]", "[CH]=C", "RadicalFeatures");
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 2};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(NoRadical) {
        FeatureVector *fv = getFeatureVector("CCC[CH4+]", "CCC", "RadicalFeatures");
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 3};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()