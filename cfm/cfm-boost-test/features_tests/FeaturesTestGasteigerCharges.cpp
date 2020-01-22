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

// Testing Gasteiger Charges Feature
BOOST_AUTO_TEST_SUITE(FeaturesTestGasteigerCharges)

// "C-C-N,C-C-N,C-X-C"
    BOOST_AUTO_TEST_CASE(GasteigerChargesTest) {
        std::vector<std::string> fnames;
        fnames.push_back("GasteigerCharges");
        FeatureCalculator *fc = new FeatureCalculator(fnames);

        romol_ptr_t ion = createMolPtr("C");
        initMolProps(ion);
        ion.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.33042488);
        RootedROMol rtd_ion(ion, ion.get()->getAtomWithIdx(0));

        romol_ptr_t nl = createMolPtr("C");
        initMolProps(nl);
        nl.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.00652530);
        RootedROMol rtd_nl(nl, nl.get()->getAtomWithIdx(0));

        FeatureVector *fv = fc->computeFeatureVector(&rtd_ion, &rtd_nl, nullptr);

        std::vector<unsigned> expected_feature_vector = {0, 9};
        // Test feature name list sizes
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fc;
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()
