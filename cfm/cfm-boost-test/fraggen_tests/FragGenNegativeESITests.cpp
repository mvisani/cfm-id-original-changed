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

struct NegativeESIFixture{
    NegativeESIFixture() {
        graph = getTestGraph("CC(=O)O", NEGATIVE_ESI_IONIZATION_MODE, true);
    }
    ~NegativeESIFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestNegativeESI, NegativeESIFixture)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 4);
    }


    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,59.01385292,"CC(=O)[O-]"},
            {1,15.02402368,"[CH3-]"},
            {2,17.00328823,"[OH-]"},
            {3,41.00328823 ,"C#C[O-]"}
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()

struct NegativeESIRingFixture{
    NegativeESIRingFixture() {
        graph = getTestGraph("C1=CN=CN=C1", NEGATIVE_ESI_IONIZATION_MODE, true);
    }
    ~NegativeESIRingFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestNegativeESIRing, NegativeESIRingFixture)

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,79.03017168,"[c-]1cncnc1"},
            {1,52.01927264,"[C-]#CN=C"},
            {2,26.00362258,"[C-]#N"},
            {3,25.00837361,"[C-]#C"},
            {4,28.01927264,"C=[N-]"},
            {5,27.02402368,"[CH-]=C"},
            {6,53.01452161,"[CH-]=NC#N"},
            {7,52.01927264,"C=[C-]C#N"},
            {8,50.00362258,"[C-]#CC#N"}
    };

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 8);
    }

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }
BOOST_AUTO_TEST_SUITE_END()

struct NegavtiveESIRingFixture{
    NegavtiveESIRingFixture() {
        graph = getTestGraph("C1=CN=CN=C1", NEGATIVE_ESI_IONIZATION_MODE, true);
    }
    ~NegavtiveESIRingFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestNegavtiveESIRing, NegavtiveESIRingFixture)

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,79.03017168,"[c-]1cncnc1"},
            {2,53.01452161,"[CH-]=NC#N"},
            {3,26.00362258,"[C-]#N"},
            {4,25.00837361,"[C-]#C"},
            {5,27.02402368,"[CH-]=C"},
            {6,28.01927264,"C=[N-]"},
            {7,50.00362258,"[C-]#CC#N"},
            {8,52.01927264,"C=[C-]C#N"},
            {10,52.01927264,"[C-]#CN=C"},
            {11,77.01452161,"[CH-]=NC#CC#N"},
    };

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 12);
    }

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

    // Ring Breaks produced by old cfmid
    // let us check if we cover those cases
    std::vector<std::tuple<int, double, std::string>> expected_ring_breaks {
            {0,79.03017168,"[c-]1cncnc1"},
            {10,52.01927264,"[C-]#CN=C"},
            {3,26.00362258,"[C-]#N"},
            {4,25.00837361,"[C-]#C"},
            {6,28.01927264,"C=[N-]"},
            {5,27.02402368,"[CH-]=C"},
            {2,53.01452161,"[CH-]=NC#N"},
            {8,52.01927264,"C=[C-]C#N"},
            {7,50.00362258,"[C-]#CC#N"},
    };

    BOOST_DATA_TEST_CASE(RingBreakFragmentTests, bdata::make(expected_ring_breaks),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
    }

    std::vector<std::tuple<int, double, std::string>> expected_intermediate_fragments {
            {1,79.03017168,"C=C=C=NC=[N-]"},
            {9,79.03017168,"CN=C=C=C=[N-]"}

    };
    BOOST_DATA_TEST_CASE(IntermediateFragmentTests, bdata::make(expected_intermediate_fragments),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
        BOOST_CHECK(graph->getFragmentAtIdx(idx)->isIntermediate());
    }

BOOST_AUTO_TEST_SUITE_END()