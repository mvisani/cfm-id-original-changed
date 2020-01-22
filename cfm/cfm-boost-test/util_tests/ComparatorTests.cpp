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
#include "Comparators.h"
#include <boost/filesystem.hpp>


void testComparator(Comparator &cmp, MolData * exact_match_moldata,
                    MolData * no_match_moldata, MolData * some_match_moldata,
                    double expected_scores[3],double max_score = 100) {

    double tol = 1e-8;
    for (int energy= 0; energy < 3; energy++){
        double score = cmp.computeScore(exact_match_moldata->getSpectrum(energy),
                                        exact_match_moldata->getPredictedSpectrum(energy));
        BOOST_CHECK_CLOSE_FRACTION(score, max_score, tol);
    }

    for (int energy= 0; energy < 3; energy++){
        double score = cmp.computeScore(some_match_moldata->getSpectrum(energy),
                                        some_match_moldata->getPredictedSpectrum(energy));
        BOOST_CHECK_CLOSE_FRACTION(score, expected_scores[energy], tol);
    }

    for (int energy= 0; energy < 3; energy++){
        double score = cmp.computeScore(no_match_moldata->getSpectrum(energy),
                                        no_match_moldata->getPredictedSpectrum(energy));
        BOOST_CHECK_CLOSE_FRACTION(score, 0.0, tol);
    }
}

struct ComparatorsTestFixture{
    ComparatorsTestFixture() {

        //Test with exact matches
        config_t cfg;
        initDefaultConfig(cfg);
        exact_match_moldata = new MolData("exact_match_moldata", "C", &cfg);
        std::string specfile = "test_data/test_spec/Test1.txt";
        exact_match_moldata->readInSpectraFromFile(specfile);
        std::string pspecfile = "test_data/test_pspec/Test1.txt";
        exact_match_moldata->readInSpectraFromFile(pspecfile, true);

        no_match_moldata = new MolData("no_match_moldata", "C", &cfg);
        specfile = "test_data/test_spec/Test2.txt";
        no_match_moldata->readInSpectraFromFile(specfile);
        pspecfile = "test_data/test_pspec/Test2.txt";
        no_match_moldata->readInSpectraFromFile(pspecfile, true);

        some_match_moldata = new MolData("some_match_moldata", "C", &cfg);
        specfile = "test_data/test_spec/Test3.txt";
        some_match_moldata->readInSpectraFromFile(specfile);
        pspecfile = "test_data/test_pspec/Test3.txt";
        some_match_moldata->readInSpectraFromFile(pspecfile, true);
    };
    ~ComparatorsTestFixture() {
        delete exact_match_moldata;
        delete no_match_moldata;
        delete some_match_moldata;
    };
    MolData * exact_match_moldata = nullptr;
    MolData * no_match_moldata = nullptr;
    MolData * some_match_moldata = nullptr;

};

BOOST_FIXTURE_TEST_SUITE(ComparatorsTests, ComparatorsTestFixture)

    BOOST_AUTO_TEST_CASE(WeightedRecallTests) {

        double tol = 1e-8;
        double expected_scores[3] = {100.0, 65.0, 50.0};
        WeightedRecall cmp(10.0, 0.01);
        testComparator(cmp, exact_match_moldata, no_match_moldata, some_match_moldata, expected_scores);
    }

    BOOST_AUTO_TEST_CASE(RecallTests) {

        double expected_scores[3] = {100.0, 66.6666666666666, 50.0};
        Recall cmp(10.0, 0.01);
        testComparator(cmp, exact_match_moldata, no_match_moldata, some_match_moldata, expected_scores);
    }

    BOOST_AUTO_TEST_CASE(PrecisionTests) {

        Precision cmp(10.0, 0.01);
        double expected_scores[3] = {50.0, 57.142857142, 50.0};

        testComparator(cmp, exact_match_moldata, no_match_moldata, some_match_moldata, expected_scores);
    }

    BOOST_AUTO_TEST_CASE(WeightedPrecisionTests) {

        WeightedPrecision cmp(10.0, 0.01);
        double expected_scores[3] = {57.142857142, 70.0, 55.5555555555555};
        testComparator(cmp, exact_match_moldata, no_match_moldata, some_match_moldata, expected_scores);
    }

    BOOST_AUTO_TEST_CASE(JaccardTests) {

        Jaccard cmp(10.0, 0.01);
        double expected_scores[3] = {0.666666667, 0.61538461538461542, 0.5};
        testComparator(cmp, exact_match_moldata, no_match_moldata, some_match_moldata, expected_scores, 1.0);
    }

BOOST_AUTO_TEST_SUITE_END()
