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


// Testing Break Atom Pair Feature
BOOST_AUTO_TEST_SUITE(FeaturesTestBreakAtomPair)

// Test CN Pair
    BOOST_AUTO_TEST_CASE(CNPair) {
        FeatureVector *fv = getFeatureVector("C", "N", "BreakAtomPair");

        // Compare Results
        std::vector<unsigned> expected_feature_vector = {0, 2};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

// Test X-C caes
// Test load from config file
    BOOST_AUTO_TEST_CASE(XCPair) {
        FeatureVector *fv = getFeatureVector("B", "C", "BreakAtomPair");
        // Compare Results
        std::vector<unsigned> expected_feature_vector = {0, 31};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()