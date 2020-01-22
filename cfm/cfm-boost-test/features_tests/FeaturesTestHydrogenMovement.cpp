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

// Testing HydrogenMovement Feature
FeatureVector *getFeatureTestHydrgeMovementFV(std::string ion_str, std::string nl_str, double h_movement) {
    std::vector<std::string> fnames;
    fnames.push_back("HydrogenMovement"); //IonRootTriples should work the same
    FeatureCalculator *fc = new FeatureCalculator(fnames);

    //Positive within range
    romol_ptr_t ion = createMolPtr(ion_str.c_str());
    initMolProps(ion);

    ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
    RootedROMol rtd_ion(ion, ion.get()->getAtomWithIdx(0));
    romol_ptr_t nl = createMolPtr(nl_str.c_str());
    initMolProps(nl);
    RootedROMol rtd_nl(nl, nl.get()->getAtomWithIdx(0));

    FeatureVector *fv = fc->computeFeatureVector(&rtd_ion, &rtd_nl, nullptr);

    delete fc;
    return fv;
}

BOOST_AUTO_TEST_SUITE(FeaturesTestHydrogenMovement)

// Postive Within Range
    BOOST_AUTO_TEST_CASE(PostiveWithinRange) {
        double h_movement = 3.00452;
        FeatureVector *fv = getFeatureTestHydrgeMovementFV("C", "C", h_movement);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 8};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

// Negative Within Range
    BOOST_AUTO_TEST_CASE(NegativeWithinRange) {
        double h_movement = -1.115;
        FeatureVector *fv = getFeatureTestHydrgeMovementFV("C", "C", h_movement);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 4};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

// Zero
    BOOST_AUTO_TEST_CASE(Zero) {
        double h_movement = 0.0;
        FeatureVector *fv = getFeatureTestHydrgeMovementFV("C", "C", h_movement);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 5};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(PositiveOutofRange) {
        double h_movement = 5.0;
        FeatureVector *fv = getFeatureTestHydrgeMovementFV("C", "C", h_movement);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 10};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

    BOOST_AUTO_TEST_CASE(NegativeOutofRange) {
        double h_movement = -5.0;
        FeatureVector *fv = getFeatureTestHydrgeMovementFV("C", "C", h_movement);

        // Check Results
        std::vector<unsigned> expected_feature_vector = {0, 10};
        BOOST_CHECK_EQUAL_COLLECTIONS(fv->getFeatureBegin(), fv->getFeatureEnd(),
                                      expected_feature_vector.begin(), expected_feature_vector.end());
        delete fv;
    }

BOOST_AUTO_TEST_SUITE_END()