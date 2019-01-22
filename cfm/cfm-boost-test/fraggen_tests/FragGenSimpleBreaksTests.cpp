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

BOOST_AUTO_TEST_SUITE(FragGenSimpleBreakTests)

    //Expected:
    //0 61.02840582 CC(O)=[OH+]
    //1 44.99710569 O=C=[OH+]
    //2 15.02292652 [CH3+]
    //3 17.03857658 [CH5+]
    //4 19.01784114 [OH3+]
    //5 43.01784114 C#C[OH2+]
    //6 25.00727645 [C+]#C
    //7 41.00219107 [C+]#CO
    //8 59.01275576 [CH+]=C(O)O

    BOOST_AUTO_TEST_CASE(FragGenTestPositiveESI) {

        std::string smiles_or_inchi = "CC(=O)O";
        int ionization_mode = POSITIVE_ESI_IONIZATION_MODE;
        FragmentGraph *graph = getTestGraph(smiles_or_inchi, ionization_mode);

        double tolerance =  1e-6;
        BOOST_CHECK_EQUAL(graph->getNumFragments(), 9);

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(0)->getMass(), 61.02840582, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(1)->getMass(), 44.99710569, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(2)->getMass(), 15.02292652, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(3)->getMass(), 17.03857658, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(4)->getMass(), 19.01784114, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(5)->getMass(), 43.01784114, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(6)->getMass(), 25.00727645, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(7)->getMass(), 41.00219107, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(8)->getMass(), 59.01275576, tolerance);
        delete graph;
    }

    BOOST_AUTO_TEST_CASE(FragGenTestNegativeESI) {

        //Expected:
        // 0 59.01385292 CC(=O)O"
        // 1 15.02402368 [CH3-]"
        // 2 17.00328823 [OH-]"
        // 3 41.00328823 C#C[O-]"
        double tolerance =  1e-6;

        std::string smiles_or_inchi = "CC(=O)O";
        int ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;

        FragmentGraph *graph = getTestGraph(smiles_or_inchi, ionization_mode);

        BOOST_CHECK_EQUAL(graph->getNumFragments(), 4);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(0)->getMass(), 59.01385292, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(1)->getMass(), 15.02402368, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(2)->getMass(), 17.00328823, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(3)->getMass(), 41.00328823, tolerance);

        delete graph;
    }

    BOOST_AUTO_TEST_CASE(FragGenTestPositiveEI) {

        //Expected:
        //0 60.02058079 CC(=[O+])O
        //1 43.98928066 O=C=[O+]
        //2 44.99710569 O=C=[OH+]
        //3 15.02292652 [CH3+]
        //4 16.03075155 [CH4+]
        //5 18.0100161 [OH2+]
        //6 19.01784114 [OH3+]
        //7 43.01784114 C#C[OH2+]
        //8 25.00727645 [C+]#C
        //9 41.00219107 [C+]#CO
        //10 42.0100161 C#C[OH+]
        //11 23.99945142 [C]#[C+]
        //12 39.99436604 [C+]#C[O]
        //13 58.00493072 [CH+]=C([O])O
        //14 59.01275576 [CH+]=C(O)O

        double tolerance =  1e-6;

        std::string smiles_or_inchi = "CC(=O)O";
        int ionization_mode = POSITIVE_EI_IONIZATION_MODE;

        FragmentGraph *graph = getTestGraph(smiles_or_inchi, ionization_mode);

        BOOST_CHECK_EQUAL(graph->getNumFragments(), 15);

        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(0)->getMass(), 60.02058079, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(1)->getMass(), 43.98928066, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(2)->getMass(), 44.99710569, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(3)->getMass(), 15.02292652, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(4)->getMass(), 16.03075155, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(5)->getMass(), 18.0100161, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(6)->getMass(), 19.01784114 , tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(7)->getMass(), 43.01784114, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(8)->getMass(), 25.00727645, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(9)->getMass(), 41.00219107, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(10)->getMass(), 42.0100161 , tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(11)->getMass(), 23.99945142 , tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(12)->getMass(), 39.99436604 , tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(13)->getMass(), 58.00493072 , tolerance);
        BOOST_CHECK_CLOSE_FRACTION(graph->getFragmentAtIdx(14)->getMass(), 59.01275576 , tolerance);

        delete graph;
    }

    BOOST_AUTO_TEST_CASE(FragGenTestPositiveEIMultibreak) {

    }
    
    BOOST_AUTO_TEST_CASE(FragGenTestPositiveESISplitCharge) {

    }

BOOST_AUTO_TEST_SUITE_END()