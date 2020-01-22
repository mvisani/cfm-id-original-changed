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

// Testing Extra Functional Groups Root Only Features
BOOST_AUTO_TEST_SUITE(FeaturesTestFunctionalGroupsRootOnly)

    BOOST_AUTO_TEST_CASE(FunctionalGroupsRootOnlyTest) {
        std::vector<std::string> fnames;
        fnames.push_back("IonFunctionalGroupRootOnlyFeatures");
        fnames.push_back("NLFunctionalGroupRootOnlyFeatures");

        std::string smiles_or_inchi("CCCCC(O)=O");
        int break_id = 2;

        FeatureVector *fv = getFeatureVector(smiles_or_inchi, fnames, break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 8, 10, 87, 324};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()