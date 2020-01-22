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

BOOST_AUTO_TEST_SUITE(FeaturesTestQuadraticFeatures)

    BOOST_AUTO_TEST_CASE(FeaturesTestQuadraticFeatures) {

        std::vector<std::string> fnames;
        fnames.push_back("BreakAtomPair");
        fnames.push_back("HydrogenMovement");
        fnames.push_back("QuadraticFeatures");
        FeatureCalculator *fc = new FeatureCalculator(fnames);

        //Simple initial vector with 3 bits set (indexes: 0,2,45)
        romol_ptr_t ion = createMolPtr("C");
        initMolProps(ion);

        RootedROMol rtd_ion(ion, ion.get()->getAtomWithIdx(0));
        double h_movement = 3.00452;
        ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
        romol_ptr_t nl = createMolPtr("N");
        initMolProps(nl);
        RootedROMol rtd_nl(nl, nl.get()->getAtomWithIdx(0));
        nl.get()->setProp("IsRingBreak", 0);

        FeatureVector *fv = fc->computeFeatureVector(&rtd_ion, &rtd_nl, nullptr);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 2, 45, 995};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());

        BOOST_CHECK(fv->getTotalLength() == 1129);
        BOOST_CHECK(fc->getNumFeatures() == 1129);

        delete fv;
        delete fc;
    }

BOOST_AUTO_TEST_SUITE_END()