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

BOOST_AUTO_TEST_SUITE(FeaturesTestLength)

    BOOST_AUTO_TEST_CASE(FeaturesTestLength) {

        //Create the feature calculator
        std::vector<std::string> names = FeatureCalculator::getValidFeatureNames();

        //Create some aribitrary input data
        FeatureCalculator full_fc( names );
        FragmentGraphGenerator fgen(&full_fc);
        std::string smiles_or_inchi("CCCCC(O)=O");

        FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
        std::vector<Break> breaks;
        node->generateBreaks(breaks, false, false);
        node->applyBreak(breaks[2], 0);	//Break Bond 2
        node->generateChildrenOfBreak(breaks[2], false);
        FragmentTreeNode *child = &(node->children[0]);
        Transition tmp_t( -1, -1, child->nl, child->ion );

        //Check all feature lengths
        std::vector<std::string>::const_iterator it = names.begin();
        for( ; it != names.end(); ++it ){

            std::vector<std::string> feature_list;
            feature_list.push_back( *it );
            FeatureCalculator fc( feature_list );

            //Compute the feature vector
            FeatureVector *fv = fc.computeFeatureVector(tmp_t.getIon(), tmp_t.getNeutralLoss(), node->ion);
            BOOST_CHECK(fv->getTotalLength() == fc.getNumFeatures());

            delete fv;
        }
    }

BOOST_AUTO_TEST_SUITE_END()