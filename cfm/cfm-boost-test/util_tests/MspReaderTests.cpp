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

struct MspReaderTestFixture{
    MspReaderTestFixture() {
        msp =  MspReader("test_data/nist2011_cutdown.msp", "NIST2011_");
        initDefaultConfig(cfg);

    };
    ~MspReaderTestFixture() {  };
    MspReader msp;
    config_t cfg;
};

BOOST_FIXTURE_TEST_SUITE(MspReaderTest, MspReaderTestFixture)

    BOOST_AUTO_TEST_CASE(SpectrumTest) {

        //Test first spectrum in msp
        double tol = 1e-5;
        MolData mol1( "NIST2011_1", "InChI=1S/H2/h1H", &cfg );
        mol1.readInSpectraFromMSP( msp );
        BOOST_CHECK_EQUAL(mol1.getNumSpectra(), 1);

        const Spectrum *spec = mol1.getSpectrum(0);
        BOOST_CHECK_EQUAL(spec->size(), 2);

        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 1.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity, 2.056903, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->mass, 2.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->intensity, 97.9430962, tol);

        MolData mol2( "NIST2011_71459", "InChI=1S/H2/h1H", &cfg );
        mol2.readInSpectraFromMSP( msp );
        BOOST_CHECK_EQUAL(mol2.getNumSpectra(), 1);
        spec = mol2.getSpectrum(0);

        BOOST_CHECK_EQUAL(spec->size(), 166);

        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 33.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity, 0.010019287, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(165)->mass, 469.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(165)->intensity, 0.0050096436, tol);

        //Test last spectrum in msp
        MolData mol3( "NIST2011_212964", "InChI=1S/H2/h1H", &cfg );
        mol3.readInSpectraFromMSP( msp );
        BOOST_CHECK_EQUAL(mol3.getNumSpectra(), 1);
        spec = mol3.getSpectrum(0);

        BOOST_CHECK_EQUAL(spec->size(), 30);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 192.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity,0.87150087, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(29)->mass, 1168.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(29)->intensity, 18.15324, tol);
    }
BOOST_AUTO_TEST_SUITE_END()


struct MultipleEnergiesMspReaderTestFixture{
    MultipleEnergiesMspReaderTestFixture() {
        msp =  MspReader("test_data/three_energies.msp", "");
        initDefaultConfig(cfg);

    };
    ~MultipleEnergiesMspReaderTestFixture() {  };
    MspReader msp;
    config_t cfg;
};

BOOST_FIXTURE_TEST_SUITE(MultipleEnergiesMspReaderTest, MultipleEnergiesMspReaderTestFixture)

    BOOST_AUTO_TEST_CASE(SpectrumTest) {

        //Test first spectrum in msp
        double tol = 1e-5;
        MolData mol1( "Test3", "InChI=1S/H2/h1H", &cfg );
        mol1.readInSpectraFromMSP( msp );
        BOOST_CHECK_EQUAL(mol1.getNumSpectra(), 3);

        const Spectrum *spec = mol1.getSpectrum(0);
        BOOST_CHECK_EQUAL(spec->size(), 2);

        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 54.07127368, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity, 10.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->mass, 76.08692374, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->intensity, 90.0, tol);

        spec = mol1.getSpectrum(1);
        BOOST_CHECK_EQUAL(spec->size(), 2);


        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 15.02292652, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity, 40.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->mass,  47.06037464, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->intensity, 60.0, tol);

        spec = mol1.getSpectrum(2);
        BOOST_CHECK_EQUAL(spec->size(), 2);

        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 15.02292652, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity,20.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->mass, 29.01342445, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->intensity, 80.0, tol);

        MolData mol2( "Test5", "InChI=1S/H2/h1H", &cfg );
        mol2.readInSpectraFromMSP( msp );
        BOOST_CHECK_EQUAL(mol2.getNumSpectra(), 3);

        spec = mol2.getSpectrum(0);
        BOOST_CHECK_EQUAL(spec->size(), 2);

        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 241.1798211, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity, 50.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->mass, 315.2166005, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(1)->intensity, 50.0, tol);

        spec = mol2.getSpectrum(1);
        BOOST_CHECK_EQUAL(spec->size(), 3);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass, 299.1853004, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity, 25.0, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(2)->mass,  315.2166005, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(2)->intensity, 25.0, tol);

        spec = mol2.getSpectrum(2);
        BOOST_CHECK_EQUAL(spec->size(), 1);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->mass,  29.03857658, tol);
        BOOST_CHECK_CLOSE_FRACTION(spec->getPeak(0)->intensity, 100.0, tol);

    }
BOOST_AUTO_TEST_SUITE_END()
