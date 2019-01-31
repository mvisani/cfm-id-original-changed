/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for Params
#
# Author: Felicity Allen, Fei Wang
# Created: November 2012, 2019
#########################################################################*/
#include <boost/test/unit_test.hpp>

#include "Config.h"
#include "FragmentTreeNode.h"
#include "MolData.h"

#pragma once


class ParamTestMol : public MolData {
public:
    ParamTestMol(config_t *cfg) : MolData("Param Test Mol", "", cfg) {
        std::vector<int> param_null_eloc;

        //Create a molecule based on what was previously in test_bn_transition.txt
        fg = new FragmentGraph();
        FeatureHelper fh;
        romol_ptr_t basic_nl = createMolPtr("C");
        fg->addToGraph(FragmentTreeNode(createMolPtr("NCCCN"), basic_nl, -1, -1, &fh, param_null_eloc), -1); //id = 0
        fg->addToGraph(FragmentTreeNode(createMolPtr("[NH4+]"), basic_nl, -1, -1, &fh, param_null_eloc),
                       0); // id = 1, 0 -> 1
        fg->addToGraph(FragmentTreeNode(createMolPtr("C=CC[NH3+]"), basic_nl, -1, -1, &fh, param_null_eloc),
                       0); //id = 2, 0 -> 2
        fg->addToGraph(FragmentTreeNode(createMolPtr("[NH4+]"), basic_nl, -1, -1, &fh, param_null_eloc), 2); // 2 -> 1
        fg->addToGraph(FragmentTreeNode(createMolPtr("C=C=[NH2+]"), basic_nl, -1, -1, &fh, param_null_eloc),
                       2);  //id = 3, 2 -> 3
        fg->addToGraph(FragmentTreeNode(createMolPtr("[CH3+]"), basic_nl, -1, -1, &fh, param_null_eloc),
                       2); //id = 4, 2->4
        fg->addToGraph(FragmentTreeNode(createMolPtr("[CH3+]"), basic_nl, -1, -1, &fh, param_null_eloc), 3); //3->4
        fg->addToGraph(FragmentTreeNode(createMolPtr("[C+]#N"), basic_nl, -1, -1, &fh, param_null_eloc),
                       3); //id = 5, 3->5
        graph_computed = true;

        //Add a simple feature vector to each transition
        unsigned int num_trans = fg->getNumTransitions();
        //fvs.resize( num_trans );
        for (unsigned int i = 0; i < num_trans; i++) {
            auto fv = new FeatureVector();
            fv->addFeature(1.0);            //Bias term
            fv->addFeatureAtIdx(0.0, 10);    //Resize to 11 (Using HydrogenMovement feature to initialise params)
            //Transition 0->2
            //Transition 2->3
            //Transition 2->4
            if (i == 1 || i == 3 || i == 4)
                fv->addFeatureAtIdx(1.0, 1);
            fg->addFeatureVectorAtIdx(i, fv);
        }
    }
};

class NNParamTestMol : public MolData {
public:
    NNParamTestMol(config_t *cfg) : MolData("NN Param Test Mol", "", cfg) {

        std::vector<int> nn_param_null_eloc;

        //Create two transitions with the feature vectors as in previous tests
        fg = new FragmentGraph();
        FeatureHelper fh;
        romol_ptr_t basic_nl = createMolPtr("C");
        fg->addToGraph(FragmentTreeNode(createMolPtr("NCCCN"), basic_nl, -1, -1, &fh, nn_param_null_eloc), -1); //id = 0
        fg->addToGraph(FragmentTreeNode(createMolPtr("[NH4+]"), basic_nl, -1, -1, &fh, nn_param_null_eloc),
                       0); // id = 1, 0 -> 1
        fg->addToGraph(FragmentTreeNode(createMolPtr("C=CC[NH3+]"), basic_nl, -1, -1, &fh, nn_param_null_eloc),
                       0); //id = 2, 0 -> 2
        graph_computed = true;

        auto fv0 = new FeatureVector();
        fv0->addFeatureAtIdx(1.0, 0);
        fv0->addFeatureAtIdx(1.0, 2);
        fv0->addFeatureAtIdx(1.0, 5);
        fg->addFeatureVectorAtIdx(0, fv0);

        auto fv1 = new FeatureVector();
        fv1 = new FeatureVector();
        fv1->addFeatureAtIdx(1.0, 0);
        fv1->addFeatureAtIdx(1.0, 3);
        fv1->addFeatureAtIdx(1.0, 5);
        fg->addFeatureVectorAtIdx(1, fv1);
    }
};
