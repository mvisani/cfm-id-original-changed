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

// Testing Extra Functional Groups Features
BOOST_AUTO_TEST_SUITE(FeaturesTestExtraFunctionalGroups)

    BOOST_AUTO_TEST_CASE(ExtraFunctionalGroupsTest) {
        std::vector<std::string> fnames;
        fnames.push_back("IonExtraFunctionalGroupFeatures");
        fnames.push_back("NLExtraFunctionalGroupFeatures");
        FeatureCalculator *fc = new FeatureCalculator(fnames);

        FragmentGraphGenerator fgen(fc);
        std::string smiles_or_inchi("C1CC1CC(C)(C)C");
        FragmentTreeNode *node = fgen.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
        std::vector<Break> breaks;
        node->generateBreaks(breaks, false, false);
        node->applyBreak(breaks[0], 0);    //Break Bond 2
        node->generateChildrenOfBreak(breaks[0], false);

        FragmentTreeNode *child = &(node->children[0]);
        Transition tmp_t(-1, -1, child->nl, child->ion);
        FeatureVector *fv = fc->computeFeatureVector(tmp_t.getIon(), tmp_t.getNeutralLoss(), nullptr);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 3, 15, 30, 35, 45};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fc;
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()