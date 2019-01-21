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
#include <boost/test/unit_test.hpp>

#include "InferenceTestsCaseClasses.h"

#include "Config.h"
#include "FragmentTreeNode.h"
#include "MolData.h"
#include "Inference.h"

#include <boost/filesystem.hpp>
#include <GraphMol/RDKitBase.h>

BOOST_AUTO_TEST_SUITE(InferenceTestSpectrumMessageWithIsotopes)

    BOOST_AUTO_TEST_CASE(InferenceTestSpectrumMessageWithIsotopes) {

        double tolerance = 0.00001;

        config_t cfg;
        initDefaultConfig(cfg);

        IsotopeSpectrumTestMol moldata(&cfg);

        Inference infer(&moldata, &cfg);
        Message msg, down_msg;
        down_msg.reset(2);
        for (int i = 0; i < 2; i++)
            down_msg.addToIdx(i, 0.0);
        infer.createSpectrumMessage(msg, 0, down_msg);

        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(0)), 0.567333, tolerance);
        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(1)), 0.432667, tolerance);
    }

BOOST_AUTO_TEST_SUITE_END()
