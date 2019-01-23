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
namespace bdata = boost::unit_test::data;

#include "FragGenTestsUtils.h"

struct PostiveEIFixture{
    PostiveEIFixture() {
        graph = getTestGraph("CC(=O)O", POSITIVE_EI_IONIZATION_MODE, true);
    }
    ~PostiveEIFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPostiveEI, PostiveEIFixture)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 15);
    }

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,60.02058079,"CC(=[O+])O"},
            {1,43.98928066,"O=C=[O+]"},
            {2,44.99710569,"O=C=[OH+]"},
            {3,15.02292652,"[CH3+]"},
            {4,16.03075155,"[CH4+]"},
            {5,18.0100161,"[OH2+]"},
            {6,19.01784114,"[OH3+]"},
            {7,43.01784114,"C#C[OH2+]"},
            {8,25.00727645,"[C+]#C"},
            {9,41.00219107,"[C+]#CO"},
            {10,42.0100161,"C#C[OH+]"},
            {11,23.99945142,"[C]#[C+]"},
            {12,39.99436604,"[C+]#C[O]"},
            {13,58.00493072,"[CH+]=C([O])O"},
            {14,59.01275576,"[CH+]=C(O)O"}
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()

struct PositiveEIMultibreakFixture{
    PositiveEIMultibreakFixture() {
        graph = getTestGraph("CCNC", POSITIVE_EI_IONIZATION_MODE, true);
    }
    ~PositiveEIMultibreakFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveEIMultibreak, PositiveEIMultibreakFixture)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 30);
    }

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,59.07349929,"CC[NH+]C"},
            {1,57.05784922,"C=C[NH+]C"},
            {2,55.04219916,"C#C[NH+]C"},
            {3,56.05002419,"C#C[NH2+]C"},
            {4,41.0265491,"[CH]=[N+]=C"},
            {5,42.03437413,"C=[N+]=C"},
            {6,15.0234751,"[CH3+]"},
            {7,17.03912516,"[CH5+]"},
            {8,16.03130013,"[CH4+]"},
            {9,31.04219916,"C[NH2+]"},
            {10,32.05002419,"C[NH3+]"},
            {11,29.0265491,"C=[NH+]"},
            {12,30.03437413,"C=[NH2+]"},
            {13,27.01089903,"C#[N+]"},
            {14,28.01872406,"C#[NH+]"},
            {15,27.0234751,"C#[CH2+]"},
            {16,26.01565006,"[CH]=[CH+]"},
            {17,29.03912516,"C=[CH3+]"},
            {18,28.03130013,"[CH2][CH2+]"},
            {19,31.05477522,"C[CH4+]"},
            {20,30.04695019,"C[CH3+]"},
            {21,42.03437413,"C#C[NH3+]"},
            {22,41.0265491,"C#C[NH2+]"},
            {23,58.06567426,"C=C[NH2+]C"},
            {24,43.04219916,"C=[N+]C"},
            {25,44.05002419,"C=[NH+]C"},
            {26,44.05002419,"C=C[NH3+]"},
            {27,18.03437413,"[NH4+]"},
            {28,43.04219916,"C=C[NH2+]"},
            {29,17.0265491,"[NH3+]"}
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()