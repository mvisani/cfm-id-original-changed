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

// Testing Root Triples Feature
// IonRootTriples should work the same
BOOST_AUTO_TEST_SUITE(FeaturesTestRootTriples)

// "C-C-N,C-C-N,C-X-C"
    BOOST_AUTO_TEST_CASE(CCNCCNCXCTriples) {
        FeatureVector *fv = getFeatureVector("C", "C(BC)C(N)N", "NLRootTriples");

        // Test feature name list sizes
        std::vector<unsigned> expected_feature_vector = {0, 6, 7, 122};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

// "N-X-X"
    BOOST_AUTO_TEST_CASE(NXXTriples) {
        FeatureVector *fv = getFeatureVector("C", "NBB", "NLRootTriples");

        // Test feature name list sizes
        std::vector<unsigned> expected_feature_vector = {0, 286};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()