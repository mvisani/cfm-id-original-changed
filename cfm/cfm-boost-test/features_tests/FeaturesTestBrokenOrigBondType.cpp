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


BOOST_AUTO_TEST_SUITE(FeaturesTestBrokenOrigBondType)

    BOOST_AUTO_TEST_CASE(IONIC) {
        std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
        int break_id = 0;
        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "BrokenOrigBondType", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 6};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(HOnly) {
        std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
        int break_id = 31;
        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "BrokenOrigBondType", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 7};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(SINGLE) {
        std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
        int break_id = 2;
        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "BrokenOrigBondType", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 1};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(DOUBLE) {
        std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
        int break_id = 5;
        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "BrokenOrigBondType", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 2};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(TRIPLE) {
        std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
        int break_id = 1;
        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "BrokenOrigBondType", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 3};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(CONJUGATED) {
        std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
        int break_id = 5;
        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "BrokenOrigBondType", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {8, 5};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(AROMATIC) {
        std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
        int break_id = 19;
        FeatureVector *fv = getFeatureVector(smiles_or_inchi, "BrokenOrigBondType", break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 4};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()
