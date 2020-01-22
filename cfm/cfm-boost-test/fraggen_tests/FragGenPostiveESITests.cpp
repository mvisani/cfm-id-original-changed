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

struct PostiveESIFixture{
    PostiveESIFixture() {
        graph = getTestGraph("CC(=O)O", POSITIVE_ESI_IONIZATION_MODE, true);
    }
    ~PostiveESIFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPostiveESI, PostiveESIFixture)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 9);
    }

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,61.02840582,"CC(O)=[OH+]"},
            {1,44.99710569,"O=C=[OH+]"},
            {2,15.02292652,"[CH3+]"},
            {3,17.03857658,"[CH5+]"},
            {4,19.01784114,"[OH3+]"},
            {5,43.01784114,"C#C[OH2+]"},
            {6,25.00727645,"[C+]#C"},
            {7,41.00219107,"[C+]#CO"},
            {8,59.01275576,"[CH+]=C(O)O"}
    };

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()

struct PositiveESIRingFixture{
    PositiveESIRingFixture() {
        graph = getTestGraph("C1=CN=CN=C1", POSITIVE_ESI_IONIZATION_MODE, true);
    }
    ~PositiveESIRingFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveESIRing, PositiveESIRingFixture)

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,81.04472458,"c1cnc[nH+]c1"},
            {2,55.02907452,"C=[NH+]C#N"},
            {3,28.01817548,"C#[NH+]"},
            {4,26.00252542,"[C+]#N"},
            {5,53.01342445,"C#[N+]C#N"},
            {6,25.00727645,"[C+]#C"},
            {7,27.02292652,"C#[CH2+]"},
            {8,29.03857658,"C=[CH3+]"},
            {9,39.02292652,"[C+]#CC"},
            {10,30.03382555,"C=[NH2+]"},
            {11,50.00252542,"[C+]#CC#N"},
            {12,52.01817548,"C#CC#[NH+]"},
            {13,54.03382555,"C=CC#[NH+]"},
            {14,65.01342445,"[C+]#CN=C=N"},
            {15,15.02292652,"[CH3+]"},
            {16,79.02907452,"[CH2+]C#CN=C=N"},
            {17,77.01342445,"[CH+]=C=C=NC#N"},
            {19,40.01817548,"[CH2+]C#N"},
            {20,65.01342445,"N#C[C+]=C=N"},
            {21,54.03382555,"C#C[NH+]=C"},
            {22,52.01817548,"C#C[N+]#C"},
            {23,79.02907452,"C=[NH+]C#CC#N"},
            {24,77.01342445,"C#[N+]C#CC#N"}
    };

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 25);
    }

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

    std::vector<std::tuple<int, double, std::string>> expected_ring_breaks {
            {0,81.04472458,"c1cnc[nH+]c1"},
            {20,65.01342445,"N#C[C+]=C=N"},
            {15,15.02292652,"[CH3+]"},
            {21,54.03382555,"C#C[NH+]=C"},
            {3,28.01817548,"C#[NH+]"},
            {4,26.00252542,"[C+]#N"},
            {6,25.00727645,"[C+]#C"},
            {7,27.02292652,"C#[CH2+]"},
            {12,52.01817548,"C#CC#[NH+]"},
            {10,30.03382555,"C=[NH2+]"},
            {19,40.01817548,"[CH2+]C#N"},
            {8,29.03857658,"C=[CH3+]"},
            {5,53.01342445,"C#[N+]C#N"},
            {2,55.02907452,"C=[NH+]C#N"},
            {13,54.03382555,"C=CC#[NH+]"},
            {12,52.01817548,"C#CC#[NH+]"},
            {11,50.00252542,"[C+]#CC#N"},
            {9,39.02292652,"[C+]#CC"},
            {14,65.01342445,"[C+]#CN=C=N"},
    };

    BOOST_DATA_TEST_CASE(RingBreakFragmentTests, bdata::make(expected_ring_breaks),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
    }

    std::vector<std::tuple<int, double, std::string>> expected_intermediate_fragments {
            {1,81.04472458,"C=C=C=[NH+]C=N"},
            {18,81.04472458,"C[NH+]=C=C=C=N"},
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

struct PositiveESISplitChargeFixture{
    PositiveESISplitChargeFixture() {
        graph = getTestGraph("CC=[N+]=[N-]", POSITIVE_ESI_IONIZATION_MODE, true);
    }
    ~PositiveESISplitChargeFixture() { delete graph; };
    FragmentGraph *graph;
};

BOOST_FIXTURE_TEST_SUITE(FragGenTestPositiveESISplitCharge, PositiveESISplitChargeFixture)

    std::vector<std::tuple<int, double, std::string>> expected_value {
            {0,57.04472458,"[CH4+]C=[N+]=[N-]"},
            {1,15.02292652,"[CH3+]"},
            {2,31.02907452,"N=[NH2+]"},
            {3,29.01342445,"N#[NH+]"},
            {4,25.00727645,"[C+]#C"},
            {5,27.02292652,"C#[CH2+]"},
            {6,29.03857658,"C=[CH3+]"},
            {7,40.01817548,"[C+]#CN"},
            {8,55.02907452,"C#C[NH+]=N"},
            {9,53.01342445,"C#C[N+]#N"},
    };

    BOOST_AUTO_TEST_CASE(NumOfFragmentTest){
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 10);
    }

    BOOST_DATA_TEST_CASE(FragmentTests, bdata::make(expected_value),
                         idx, expected_mass, expected_ion_smiles)
    {
        double tolerance =  1e-6;

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(FragGenTestPositiveESIMaxElectronMovement)

    BOOST_AUTO_TEST_CASE(NumOfFragmentTestAtTopLeveL){
        FragmentGraphGenerator fgen(0);
        std::string smiles_or_inchi("C#CC#CC#CCCCCCCOCCCCCCCN");
        FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );

        //Break just the one bond
        std::vector<Break> breaks;
        node->generateBreaks(breaks, false);
        node->applyBreak(breaks[12], 0);	//Break Bond 11 (after the O)
        node->generateChildrenOfBreak(breaks[12]);

        //Creat a simple graph for just these breaks
        config_t cfg; initDefaultConfig(cfg);
        cfg.include_isotopes = false; cfg.include_h_losses = true;
        FragmentGraph fg(&cfg);
        fg.addToGraph( *node, -1 );
        std::vector<FragmentTreeNode>::iterator itt = node->children.begin();
        for( ; itt != node->children.end(); ++itt ){
            fg.addToGraph( *itt, 0 );
        }

            //Now break one of the children (to check the second level works)
            FragmentTreeNode *child = &node->children[6];
            std::vector<Break> child_breaks;
            child->generateBreaks(child_breaks, false);
            child->applyBreak(breaks[5], 0);	//Break Bond 5
            child->generateChildrenOfBreak(breaks[5]);
            itt = child->children.begin();
            for( ; itt != child->children.end(); ++itt ){
                fg.addToGraph( *itt, 0 );
            }
        BOOST_CHECK_EQUAL(fg.getNumFragments(), 20);
        delete node;
    }

BOOST_AUTO_TEST_SUITE_END()
