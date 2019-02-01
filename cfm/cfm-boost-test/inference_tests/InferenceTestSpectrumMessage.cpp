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

#include "InferenceTestMols.h"

#include "Config.h"
#include "FragmentTreeNode.h"
#include "MolData.h"
#include "Inference.h"

#include <boost/filesystem.hpp>
#include <GraphMol/RDKitBase.h>

BOOST_AUTO_TEST_SUITE(InferenceTests)

    BOOST_AUTO_TEST_CASE(InferenceTestSpectrumMessage) {

        double tolerance = 0.00001;

        config_t cfg;
        initDefaultConfig(cfg);

        SpectrumTestMol moldata(&cfg);

        Inference infer(&moldata, &cfg);
        Message msg, down_msg;
        down_msg.reset(6);
        for (int i = 0; i < 6; i++)
            down_msg.addToIdx(i, 0.0);
        infer.createSpectrumMessage(msg, 0, down_msg);

        BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(0)), 0.1, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(1)), 0.15, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(2)), 0.15, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(3)), 0.2, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(4)), 0.2, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(5)), 0.2, tolerance);
    }

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

            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(0)), 0.1, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(1)), 0.15, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(2)), 0.15, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(3)), 0.2, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(4)), 0.2, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(5)), 0.2, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(5)), 0.0, tolerance);

    }
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

            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(0)), 0.567333, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(1)), 0.432667, tolerance);
    }

    BOOST_AUTO_TEST_CASE(InferenceTestSpectrumMessageWithIsotopesAndNoisePeak) {

            double tolerance = 0.00001;

            config_t cfg;
            initDefaultConfig(cfg);
            cfg.abs_mass_tol = 0.1;
            cfg.ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;
            cfg.fg_depth = 2;
            cfg.include_isotopes = true; cfg.isotope_thresh = 0.01;

            IsotopeSpectrumTestMolNoisePeak moldata(&cfg);

            Inference infer( &moldata, &cfg);
            Message msg, down_msg;
            down_msg.reset(3); for(int i=0; i<3; i++ ) down_msg.addToIdx(i, 0.0);
            infer.createSpectrumMessage( msg, 0, down_msg );


            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(0)), 0.567333, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(1)), 0.432667, tolerance);
            BOOST_CHECK_CLOSE_FRACTION(std::exp(msg.getIdx(2)), 0.0, tolerance);
    }

BOOST_AUTO_TEST_SUITE_END()