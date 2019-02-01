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

// Testing Root Pairs Feature
// NLRootPairs should work the same
BOOST_AUTO_TEST_SUITE(FeaturesTestRootPairs)

// "C-N,C-N,C-X"
    BOOST_AUTO_TEST_CASE(CNCNCXPairs) {
        FeatureVector *fv = getFeatureVector("C(N)(B)N", "C", "IonRootPairs");
        // Compare Results
        std::vector<unsigned> expected_feature_vector = {0, 6, 7, 22};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

// "X-X,X-N"
    BOOST_AUTO_TEST_CASE(XXXNPairs) {
        FeatureVector *fv = getFeatureVector("B(B)N", "C", "IonRootPairs");
        // Compare Results
        std::vector<unsigned> expected_feature_vector = {0, 126, 142};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

// EmptyPairs
    BOOST_AUTO_TEST_CASE(EmptyPairs) {
        FeatureVector *fv = getFeatureVector("C", "C", "IonRootPairs");
        // Compare Results
        std::vector<unsigned> expected_feature_vector = {0, 1};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()