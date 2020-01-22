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
#include "FragGenTestsUtils.h"

FragmentGraph *getTestGraph(std::string smiles_or_inchi,
                            int ionization_mode,
                            bool include_h_loss,
                            bool allow_frag_detours,
                            int graph_depth) {
    //Run the fragmentation procedure
    FragmentGraphGenerator gg;
    FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, ionization_mode);

    config_t cfg;
    cfg.ionization_mode = ionization_mode;
    initDefaultConfig(cfg);
    cfg.include_h_losses = include_h_loss;
    cfg.allow_frag_detours = allow_frag_detours;
    FragmentGraph *graph = gg.createNewGraph(&cfg);

    gg.compute(*startNode, graph_depth, -1, 2);
    delete startNode;
    return graph;
}