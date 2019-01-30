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
#include <boost/test/data/test_case.hpp>

#include "MolData.h"
#include <boost/filesystem.hpp>

struct SpectrumQuantiseTestFixture{
    SpectrumQuantiseTestFixture() {

        config_t cfg; initDefaultConfig(cfg);
        moldata = new MolData ("Test4","C", &cfg);
        std::string specfile = "test_data/test_pspec/Test3.txt";

        moldata->readInSpectraFromFile(specfile, true);
    };
    ~SpectrumQuantiseTestFixture() { delete moldata; };
    MolData * moldata = nullptr;
};

BOOST_FIXTURE_TEST_SUITE(SpectrumQuantiseTest, SpectrumQuantiseTestFixture)

    BOOST_AUTO_TEST_CASE(ThreeDecimal) {
        double tol = 1e-10;
        moldata->quantisePredictedSpectra(3);
        const Spectrum *spec = moldata->getPredictedSpectrum(1);

        double exp_masses[7] = {10.0, 50.99, 60.0, 76.008, 76.31, 76.431, 76.432};
        double exp_intensities[7] = {20.0, 20.0, 35.0, 5.0, 5.0, 5.0, 10.0};

        BOOST_CHECK_EQUAL(spec->size(), 7);

        for (int i = 0; i < 7; i++) {
            BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(i)->mass, exp_masses[i], tol);
            BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(i)->intensity, exp_intensities[i], tol);
        }
    }

    BOOST_AUTO_TEST_CASE(TwoDecimal) {
        double tol = 1e-10;

        moldata->quantisePredictedSpectra(2);
        const Spectrum *spec = moldata->getPredictedSpectrum(1);

        double exp_masses[6] = {10.0, 50.99, 60.0, 76.01, 76.31, 76.43};
        double exp_intensities[6] = {20.0, 20.0, 35.0, 5.0, 5.0, 15.0};

        BOOST_CHECK_EQUAL(spec->size(), 6);

        for (int i = 0; i < 6; i++) {
            BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(i)->mass, exp_masses[i], tol);
            BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(i)->intensity, exp_intensities[i], tol);
        }

    }

    BOOST_AUTO_TEST_CASE(ZeroDecimal) {
        double tol = 1e-10;

        moldata->quantisePredictedSpectra(0);
        const Spectrum *spec = moldata->getPredictedSpectrum(1);

        double exp_masses[4] = {10.0, 51.0, 60.0, 76.0};
        double exp_intensities[4] = {20.0, 20.0, 35.0, 25.0};


        BOOST_CHECK_EQUAL(spec->size(), 4);

        for (int i = 0; i < 4; i++) {
            BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(i)->mass, exp_masses[i], tol);
            BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(i)->intensity, exp_intensities[i], tol);
        }
    }

BOOST_AUTO_TEST_SUITE_END()


/* Disabled Since missing this test file
struct SpectrumCleanTestFixture{
    SpectrumCleanTestFixture() {

        config_t cfg; initDefaultConfig(cfg);
        moldata = new MolData ("Test4,"C", &cfg);
        std::string specfile = "test_data/test_pspec/Test4.txt";

        moldata->readInSpectraFromFile(specfile, true);
        moldata->cleanSpectra(0.1, 10.0);
    };
    ~SpectrumCleanTestFixture() { delete moldata; };
    MolData * moldata = nullptr;
};

BOOST_FIXTURE_TEST_SUITE(SpectrumCleanTest, SpectrumCleanTestFixture)

    BOOST_AUTO_TEST_CASE(cleanTest) {

        for( int energy = 0; energy < 3; ++energy){
            const Spectrum *spec = moldata->getSpectrum(energy);
            int count = 0;
            for(auto  it = spec->begin(); it != spec->end(); ++it )
                if( fabs(it->mass - 113.07) < 0.01 && it->intensity > 10.0 )
                    count++;
            BOOST_TEST(count == 1);
        }
    }

BOOST_AUTO_TEST_SUITE_END()
 */