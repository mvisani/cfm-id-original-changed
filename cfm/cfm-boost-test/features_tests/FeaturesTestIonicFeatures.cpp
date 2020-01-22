/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code
#
# Author: Felicity Allen, Fei Wang
# Created: November 2012, 2019
#########################################################################*/
#include "FeaturesTestsUtil.h"

FeatureVector *getFeaturesTestIonicFeaturesFV(int break_idx, int ionic_idx, bool rebreak_child = false) {
    std::vector<std::string> fnames;
    fnames.push_back("IonicFeatures");
    FeatureCalculator *fc = new FeatureCalculator( fnames );

    FragmentGraphGenerator fgen(fc);
    std::string smiles_or_inchi("[Na+].[Cl-].CC=CC");
    FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
    std::vector<Break> breaks;
    node->generateBreaks(breaks, true);

    Break *brk = &breaks[break_idx];

    node->applyBreak(*brk, ionic_idx );	//Break specified bond
    node->generateChildrenOfBreak(*brk);
    FragmentTreeNode *child = &(node->children[0]);

    if( rebreak_child ){
        std::vector<Break> child_breaks;
        child->generateBreaks(child_breaks, true);
        child->applyBreak(child_breaks[0], 0 );
        child->generateChildrenOfBreak(child_breaks[0]);
        child = &(child->children[0]);
    }

    Transition tmp_t( -1, -1, child->nl, child->ion );
    FeatureVector *fv = fc->computeFeatureVector(tmp_t.getIon(), tmp_t.getNeutralLoss(), nullptr);

    return fv;
}

BOOST_AUTO_TEST_SUITE(FeaturesTestIonicFeatures)

    BOOST_AUTO_TEST_CASE(NaIonic) {
        int break_id = 0;
        int ion_idx = 0;

        FeatureVector *fv = getFeaturesTestIonicFeaturesFV(break_id, ion_idx);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 2, 3};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    //Tests case with Na+ and Cl- on NL side, Ion On CC=CC
    BOOST_AUTO_TEST_CASE(CCCCIonic) {
        int break_id = 0;
        int ion_idx = 1;

        FeatureVector *fv = getFeaturesTestIonicFeaturesFV(break_id, ion_idx);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 1, 3};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    //Tests case with Na+ and Cl- on Ion side
    BOOST_AUTO_TEST_CASE(NaClIon) {
        int break_id = 2;
        int ion_idx = 0;

        FeatureVector *fv = getFeaturesTestIonicFeaturesFV(break_id, ion_idx);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 2, 4};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    //Tests case with no ionic fragments (need to re-break a child from above)
    BOOST_AUTO_TEST_CASE(NoIonic) {
        int break_id = 0;
        int ion_idx = 1;

        FeatureVector *fv = getFeaturesTestIonicFeaturesFV(break_id, ion_idx, true);
        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 5};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }
BOOST_AUTO_TEST_SUITE_END()