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

#include "MolData.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/filesystem.hpp>


FragmentGraph *getTestGraph(std::string smiles_or_inchi, int ionization_mode) {
    //Run the fragmentation procedure
    FragmentGraphGenerator gg;
    FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, ionization_mode);

    config_t cfg;
    cfg.ionization_mode = ionization_mode;
    initDefaultConfig(cfg);
    cfg.include_h_losses = true;
    FragmentGraph *graph = gg.createNewGraph(&cfg);

    gg.compute(*startNode, 2, -1, 2);
    delete startNode;
    return graph;
}

void checkFragments(FragmentGraph * graph,
                    std::vector<std::tuple<int, double, std::string>> & expected,
                    double tolerance = 1e-5){

    for(const auto & row : expected){
        int idx = std::get<0>(row);
        double expected_mass = std::get<1>(row);
        std::string expected_ion_smiles = std::get<2>(row);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }
}

BOOST_AUTO_TEST_SUITE(FragGenSimpleBreakTests)

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

    BOOST_DATA_TEST_CASE(FragGenTestPositiveESI, bdata::make(expected_value),
            idx, expected_mass, expected_ion_smiles)
    {
        std::string smiles_or_inchi = "CC(=O)O";
        int ionization_mode = POSITIVE_ESI_IONIZATION_MODE;
        FragmentGraph *graph = getTestGraph(smiles_or_inchi, ionization_mode);

        double tolerance =  1e-6;
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 9);

        std::cout << idx << expected_mass << expected_ion_smiles << std::endl;
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(idx)->getMass(), expected_mass, tolerance);
        BOOST_CHECK_EQUAL(*(graph->getFragmentAtIdx(idx)->getIonSmiles()), expected_ion_smiles);
    }

    BOOST_AUTO_TEST_CASE(FragGenTestNegativeESI) {

        double tolerance =  1e-5;

        std::string smiles_or_inchi = "CC(=O)O";
        int ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;

        FragmentGraph *graph = getTestGraph(smiles_or_inchi, ionization_mode);

        BOOST_CHECK_EQUAL(graph->getNumFragments(), 4);
        std::vector<std::tuple<int, double, std::string>> expected;

        // Expected
        expected.push_back(std::make_tuple(0,59.01385292,"CC(=O)[O-]"));
        expected.push_back(std::make_tuple(1,15.02402368,"[CH3-]"));
        expected.push_back(std::make_tuple(2,17.00328823,"[OH-]"));
        expected.push_back(std::make_tuple(3,41.00328823 ,"C#C[O-]"));

        checkFragments(graph, expected, tolerance);

        delete graph;
    }

    BOOST_AUTO_TEST_CASE(FragGenTestPositiveEI) {


        double tolerance =  1e-6;

        std::string smiles_or_inchi = "CC(=O)O";
        int ionization_mode = POSITIVE_EI_IONIZATION_MODE;

        FragmentGraph *graph = getTestGraph(smiles_or_inchi, ionization_mode);

        BOOST_CHECK_EQUAL(graph->getNumFragments(), 15);

        std::vector<std::tuple<int, double, std::string>> expected;
        expected.push_back(std::make_tuple(0,60.02058079,"CC(=[O+])O"));
        expected.push_back(std::make_tuple(1,43.98928066,"O=C=[O+]"));
        expected.push_back(std::make_tuple(2,44.99710569,"O=C=[OH+]"));
        expected.push_back(std::make_tuple(3,15.02292652,"[CH3+]"));
        expected.push_back(std::make_tuple(4,16.03075155,"[CH4+]"));
        expected.push_back(std::make_tuple(5,18.0100161,"[OH2+]"));
        expected.push_back(std::make_tuple(6,19.01784114,"[OH3+]"));
        expected.push_back(std::make_tuple(7,43.01784114,"C#C[OH2+]"));
        expected.push_back(std::make_tuple(8,25.00727645,"[C+]#C"));
        expected.push_back(std::make_tuple(9,41.00219107,"[C+]#CO"));
        expected.push_back(std::make_tuple(10,42.0100161,"C#C[OH+]"));
        expected.push_back(std::make_tuple(11,23.99945142,"[C]#[C+]"));
        expected.push_back(std::make_tuple(12,39.99436604,"[C+]#C[O]"));
        expected.push_back(std::make_tuple(13,58.00493072,"[CH+]=C([O])O"));
        expected.push_back(std::make_tuple(14,59.01275576,"[CH+]=C(O)O"));
        checkFragments(graph, expected, tolerance);

        delete graph;
    }

    BOOST_AUTO_TEST_CASE(FragGenTestPositiveEIMultibreak) {

    }
    
    BOOST_AUTO_TEST_CASE(FragGenTestPositiveESISplitCharge) {

    }

BOOST_AUTO_TEST_SUITE_END()