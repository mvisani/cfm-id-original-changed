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

FeatureVector *getFeaturesTestRootAtomFV(int break_id) {
    std::vector<std::string> fnames;
    fnames.push_back("IonRootAtom");
    fnames.push_back("NLRootAtom");

    std::string smiles_or_inchi("[Br]C([Cl])(F)N(I)OPS[Se][Si]([Na])CCCCCc1ccccc1");

    FeatureVector *fv = getFeatureVector(smiles_or_inchi, fnames, break_id, true);
    return fv;
}

BOOST_AUTO_TEST_SUITE(FeaturesTestRootAtom)
    BOOST_AUTO_TEST_CASE(CH) {
        int break_id = 23;

        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 2, 11 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(CBr) {
        int break_id = 0;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 2, 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(CCl) {
        int break_id = 1;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 3, 15};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(CF) {
        int break_id = 2;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 3 + 1, 1 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(CN) {
        int break_id = 3;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 5 + 1, 1 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(NI) {
        int break_id = 4;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 4 + 1, 5 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(NO) {
        int break_id = 5;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 6 + 1, 5 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(OP) {
        int break_id = 6;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 7 + 1, 6 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(PS) {
        int break_id = 7;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 8 + 1, 7 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(SSe) {
        int break_id = 8;

        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 9 + 1, 8 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(SeSi) {
        int break_id = 9;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 10 + 1, 9 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(SiNa) {
        int break_id = 10;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 12 + 1, 10 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
    BOOST_AUTO_TEST_CASE(cccc) {
        int break_id = 19;
        FeatureVector *fv = getFeaturesTestRootAtomFV(break_id);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 1 + 1 , 1 + 14};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
BOOST_AUTO_TEST_SUITE_END()