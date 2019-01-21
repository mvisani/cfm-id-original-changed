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

// Testing case for InitialisationOfFeatureClass
BOOST_AUTO_TEST_SUITE(InitialisationOfFeatureClass)
    // Test load from config file
    BOOST_AUTO_TEST_CASE(LoadFromConfigFile) {
        std::string config_filename = "./test_data/valid_feature_config.txt";
        FeatureCalculator *fc = new FeatureCalculator(config_filename);
        std::vector<std::string> fnames = fc->getFeatureNames();

        std::vector<std::string> expected_fnames = {"BreakAtomPair", "GasteigerCharges"};
        // Test feature name list sizes
        BOOST_CHECK_EQUAL_COLLECTIONS(fnames.begin(), fnames.end(),
                                      expected_fnames.begin(), expected_fnames.end());

        // one for bias
        BOOST_TEST(fc->getNumFeatures() == (37 + 36 + 1));
        delete fc;
    }

    // Test load from a not existing config file
    BOOST_AUTO_TEST_CASE(LoadFromValidConfigFile) {
        std::string config_filename = "tests/test_data/invalid_feature_config.txt";
        FeatureCalculator *fc = nullptr;
        BOOST_CHECK_THROW(fc = new FeatureCalculator(config_filename), InvalidConfigException);
        if (fc)
            delete fc;
    }

BOOST_AUTO_TEST_SUITE_END()