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

BOOST_AUTO_TEST_SUITE(FeaturesTestNeighbourMMFFAtomType)

    BOOST_AUTO_TEST_CASE(NeighbourMMFFAtomTypeTest) {
        std::vector<std::string> fnames;
        fnames.push_back("IonNeighbourMMFFAtomType");
        fnames.push_back("NLNeighbourMMFFAtomType");

        std::string smiles_or_inchi("CCCCC(O)=O");
        int break_id = 4;

        FeatureVector *fv = getFeatureVector(smiles_or_inchi, fnames, break_id, false);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 101, 102, 108};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()