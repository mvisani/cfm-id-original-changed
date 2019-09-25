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
// Testing Functional Groups Features
BOOST_AUTO_TEST_SUITE(FeaturesTestFunctionalGroups)

    BOOST_AUTO_TEST_CASE(FunctionalGroupsTest) {
        std::vector<std::string> fnames;
        fnames.push_back("IonFunctionalGroupFeatures");
        fnames.push_back("NLFunctionalGroupFeatures");
        FeatureCalculator *fc = new FeatureCalculator(fnames);

        FragmentGraphGenerator fgen(fc);
        std::string smiles_or_inchi("CCCCC(O)=O");
        FragmentTreeNode *node = fgen.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
        std::vector<Break> breaks;
        node->generateBreaks(breaks, false);
        node->applyBreak(breaks[2], 0);    //Break Bond 2
        node->generateChildrenOfBreak(breaks[2]);

        FragmentTreeNode *child = &(node->children[0]);
        Transition tmp_t(-1, -1, child->nl, child->ion);
        FeatureVector *fv = fc->computeFeatureVector(tmp_t.getIon(), tmp_t.getNeutralLoss(), nullptr);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 8, 10, 87, 163, 173, 246, 270, 486};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fc;
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()

