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

BOOST_AUTO_TEST_SUITE(InferenceTestSpectrumNoisePeakMessage)

    BOOST_AUTO_TEST_CASE(InferenceTestSpectrumNoisePeakMessage) {

        double tolerance = 0.00001;

        config_t cfg;
        initDefaultConfig(cfg);

        SpectrumTestMolNoisePeak moldata(&cfg);

        Inference infer(&moldata, &cfg);
        Message msg, down_msg;
        down_msg.reset(7);
        for (int i = 0; i < 7; i++)
            down_msg.addToIdx(i, 0.0);
        infer.createSpectrumMessage(msg, 0, down_msg);

        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(0)), 0.1, tolerance);
        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(1)), 0.15, tolerance);
        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(2)), 0.15, tolerance);
        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(3)), 0.2, tolerance);
        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(4)), 0.2, tolerance);
        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(5)), 0.2, tolerance);
        BOOST_CHECK_CLOSE(std::exp(msg.getIdx(5)), 0.0, tolerance);

    }

BOOST_AUTO_TEST_SUITE_END()