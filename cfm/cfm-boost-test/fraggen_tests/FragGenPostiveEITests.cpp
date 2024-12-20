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
            {0,59.07295071,"CC[NH+]C"},
            {1,43.04165058,"C=[N+]C"},
            {2,27.01035045,"C#[N+]"},
            {3,28.01817548,"C#[NH+]"},
            {4,15.02292652,"[CH3+]"},
            {5,17.03857658,"[CH5+]"},
            {6,16.03075155,"[CH4+]"},
            {7,41.02600052,"[CH]=[N+]=C"},
            {8,42.03382555,"C=[N+]=C"},
            {9,44.04947561,"C=[NH+]C"},
            {10,31.04165058,"C[NH2+]"},
            {11,29.02600052,"C=[NH+]"},
            {12,30.03382555,"C=[NH2+]"},
            {13,32.04947561,"C[NH3+]"},
            {14,29.03857658,"C=[CH3+]"},
            {15,27.02292652,"C#[CH2+]"},
            {16,28.03075155,"[CH2][CH2+]"},
            {17,26.01510148,"[CH]=[CH+]"},
            {18,31.05422664,"C[CH4+]"},
            {19,30.04640161,"C[CH3+]"},
            {20,44.04947561,"C=C[NH3+]"},
            {21,18.03382555,"[NH4+]"},
            {22,42.03382555,"C#C[NH3+]"},
            {23,43.04165058,"C=C[NH2+]"},
            {24,17.02600052,"[NH3+]"},
            {25,41.02600052,"C#C[NH2+]"},
            {26,57.05730064,"C=C[NH+]C"},
            {27,55.04165058,"C#C[NH+]C"},
            {28,56.04947561,"C#C[NH2+]C"},
            {29,58.06512568,"C=C[NH2+]C"},
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()


struct PositiveEITripleFixture{
    PositiveEITripleFixture() {
        graph = getTestGraph("CC#CC", POSITIVE_EI_IONIZATION_MODE);
    }
    ~PositiveEITripleFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveEITriple, PositiveEITripleFixture)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 14);
    }

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,54.04640161,"C[C]=[C+]C"},
            {1,38.01510148,"[C+]#C[CH2]"},
            {2,39.02292652,"[C+]#CC"},
            {3,15.02292652,"[CH3+]"},
            {4,28.03075155,"[CH2][CH2+]"},
            {5,26.01510148,"[CH]=[CH+]"},
            {6,27.02292652,"C#[CH2+]"},
            {7,29.03857658,"C=[CH3+]"},
            {8,23.99945142,"[C]#[C+]"},
            {9,25.00727645,"[C+]#C"},
            {10,52.03075155,"[CH+]=[C]C=C"},
            {11,50.01510148,"C#C[C]=[CH+]"},
            {12,51.02292652,"C#CC#[CH2+]"},
            {13,53.03857658,"C=CC#[CH2+]"},
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()

/*struct PositiveEIOxygenAromaticFixture{
    PositiveEIOxygenAromaticFixture() {
        bool include_h_losses = true;
        bool allow_detour = false;
        int graph_depth = 3;
        graph = getTestGraph("C1=CC2OC=CC12", POSITIVE_EI_IONIZATION_MODE,
                             include_h_losses, allow_detour, graph_depth);
    }
    ~PositiveEIOxygenAromaticFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveEIOxygenAromatic, PositiveEIOxygenAromaticFixture)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 43);
    }

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,94.04131623,"C1=CC2[O+]C=CC12"},
            {1,26.01510148,"[CH]=[CH+]"},
            {2,23.99945142,"[C]#[C+]"},
            {3,25.00727645,"[C+]#C"},
            {4,27.02292652,"C#[CH2+]"},
            {5,69.0334912,"C1=C[OH+]C=C1"},
            {6,19.01784114,"[OH3+]"},
            {7,51.02292652,"C#CC#[CH2+]"},
            {8,29.00219107,"C#[O+]"},
            {9,39.02292652,"C#C[CH2+]"},
            {10,41.03857658,"C#C[CH4+]"},
            {11,43.01784114,"C=C=[OH+]"},
            {12,41.00219107,"[CH+]=C=O"},
            {13,53.00219107,"[CH+]=C=C=O"},
            {14,15.02292652,"[CH3+]"},
            {15,17.03857658,"[CH5+]"},
            {16,68.02566617,"C1=C[O+]C=C1"},
            {17,18.0100161,"[OH2+]"},
            {18,50.01510148,"C#C[C]=[CH+]"},
            {19,27.99436604,"[C]#[O+]"},
            {20,38.01510148,"[C]#C[CH2+]"},
            {21,40.03075155,"[CH+]=[C]C"},
            {22,42.0100161,"C=C=[O+]"},
            {23,16.03075155,"[CH4+]"},
            {24,51.99436604,"[C]#CC#[O+]"},
            {25,53.03857658,"C1=C[CH2+]=C1"},
            {26,52.03075155,"[CH]1C=C[CH+]1"},
            {27,65.03857658,"[CH+]=C1C=CC1"},
            {28,49.00727645,"[C+]#CC#C"},
            {29,64.03075155,"[CH+]=C1[C]=CC1"},
            {30,47.99945142,"[C]#CC#[C+]"},
            {31,66.04640161,"C=C1[CH+][CH]C1"},
            {32,77.03857658,"[CH2+]#CC1=CC=C1"},
            {33,75.02292652,"[C+]#CC1=CC=C1"},
            {34,76.03075155,"C#CC1=C[CH][CH+]1"},
            {35,74.01510148,"[C+]#CC1=CC=[C]1"},
            {36,68.02566617,"[O+]=C1C=CC1"},
            {37,39.99436604,"[C+]#C[O]"},
            {38,69.0334912,"[OH+]=C1C=CC1"},
            {39,78.0100161,"[CH+]=C1[C]=CC1=O"},
            {40,79.01784114,"[CH+]=C1C=CC1=O"},
            {41,92.02566617,"C1=CC2=C1C=C[O+]2"},
            {42,93.0334912,"C1=CC2=C1C=C[OH+]2"}
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()*/

struct PositiveEIRingFixture{
    PositiveEIRingFixture() {
        graph = getTestGraph("C1=CN=CN=C1", POSITIVE_EI_IONIZATION_MODE, true);
    }
    ~PositiveEIRingFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveEIRing, PositiveEIRingFixture)

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,80.03689955,"C1=C[N+]=CN=C1"},
            {2,54.02124948,"C=[N+]C#N"},
            {3,27.01035045,"C#[N+]"},
            {4,28.01817548,"C#[NH+]"},
            {5,26.00252542,"[C+]#N"},
            {6,52.00559942,"[C]#[N+]C#N"},
            {7,53.01342445,"C#[N+]C#N"},
            {8,55.02907452,"C=[NH+]C#N"},
            {9,25.00727645,"[C+]#C"},
            {10,23.99945142,"[C]#[C+]"},
            {11,27.02292652,"C#[CH2+]"},
            {12,26.01510148,"[CH]=[CH+]"},
            {13,28.03075155,"[CH2][CH2+]"},
            {14,39.02292652,"[C+]#CC"},
            {15,38.01510148,"[C+]#C[CH2]"},
            {16,29.02600052,"C=[NH+]"},
            {17,30.03382555,"C=[NH2+]"},
            {18,50.00252542,"[C+]#CC#N"},
            {19,52.01817548,"C#CC#[NH+]"},
            {20,51.01035045,"C#CC#[N+]"},
            {21,54.03382555,"C=CC#[NH+]"},
            {22,53.02600052,"C=CC#[N+]"},
            {23,64.00559942,"[C+]#CN=C=[N]"},
            {24,65.01342445,"[C+]#CN=C=N"},
            {25,15.02292652,"[CH3+]"},
            {26,78.02124948,"[CH2+]C#CN=C=[N]"},
            {27,77.01342445,"[CH+]=C=C=NC#N"},
            {28,79.02907452,"[CH2+]C#CN=C=N"},
            {30,40.01817548,"[CH2+]C#N"},
            {31,65.01342445,"N#C[C+]=C=N"},
            {32,64.00559942,"[N]=C=[C+]C#N"},
            {33,53.02600052,"C#C[N+]=C"},
            {34,51.01035045,"[C]#C[N+]#C"},
            {35,52.01817548,"C#C[N+]#C"},
            {36,54.03382555,"C#C[NH+]=C"},
            {37,78.02124948,"C=[N+]C#CC#N"},
            {38,76.00559942,"[C]#[N+]C#CC#N"},
            {39,77.01342445,"C#[N+]C#CC#N"},
            {40,79.02907452,"C=[NH+]C#CC#N"}
    };

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 41);
    }

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

    std::vector<std::tuple<int, double, std::string>> expected_intermediate_fragments {
            {29,80.03689955,"C[N+]=C=C=C=N"},
            {1,80.03689955,"C=C=C=[N+]C=N"}
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

struct PositiveEIAlkaneFixture{
    PositiveEIAlkaneFixture() {
        graph = getTestGraph("CCC", POSITIVE_EI_IONIZATION_MODE, true);
    }
    ~PositiveEIAlkaneFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveEIAlkane, PositiveEIAlkaneFixture)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 12);
    }

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,44.06205168,"CC[CH3+]"},
            {1,28.03075155,"[CH2][CH2+]"},
            {2,26.01510148,"[CH]=[CH+]"},
            {3,27.02292652,"C#[CH2+]"},
            {4,29.03857658,"C=[CH3+]"},
            {5,15.02292652,"[CH3+]"},
            {6,17.03857658,"[CH5+]"},
            {7,16.03075155,"[CH4+]"},
            {8,42.04640161,"[CH2+][CH]C"},
            {9,40.03075155,"[CH+]=[C]C"},
            {10,41.03857658,"[CH2+]#CC"},
            {11,43.05422664,"CC=[CH3+]"}
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()

struct PositiveEISplitChargeFixture{
    PositiveEISplitChargeFixture() {
        graph = getTestGraph("CC=[N+]=[N-]", POSITIVE_EI_IONIZATION_MODE, true);
    }
    ~PositiveEISplitChargeFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveEISplitCharge, PositiveEISplitChargeFixture)

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,56.03689955,"[CH3+]C=[N+]=[N-]"},
            {1,15.02292652,"[CH3+]"},
            {2,30.02124948,"N=[NH+]"},
            {3,28.00559942,"N#[N+]"},
            {4,29.01342445,"N#[NH+]"},
            {5,31.02907452,"N=[NH2+]"},
            {6,25.00727645,"[C+]#C"},
            {7,23.99945142,"[C]#[C+]"},
            {8,27.02292652,"C#[CH2+]"},
            {9,26.01510148,"[CH]=[CH+]"},
            {10,28.03075155,"[CH2][CH2+]"},
            {11,40.01817548,"[C+]#CN"},
            {12,39.01035045,"[C+]#C[NH]"},
            {13,54.02124948,"C#C[N+]=N"},
            {14,52.00559942,"[C]#C[N+]#N"},
            {15,53.01342445,"C#C[N+]#N"},
            {16,55.02907452,"C#C[NH+]=N"},
    };

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 17);
    }

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()