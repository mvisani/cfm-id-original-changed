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

BOOST_AUTO_TEST_SUITE(FeaturesTestRootMMFFAtomType)

    BOOST_AUTO_TEST_CASE(RootMMFFAtomTypeTest) {
        std::vector<std::string> fnames;
        fnames.push_back("IonRootMMFFAtomType");
        fnames.push_back("NLRootMMFFAtomType");

        std::string smiles_or_inchi("CCCCC(O)=O");
        int break_id = 4;

        FeatureVector *fv = getFeatureVector(smiles_or_inchi, fnames, break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 6, 103};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()