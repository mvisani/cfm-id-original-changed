/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NLRootMatrixFPN10.cpp
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see
param.cpp.
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#################################################################y########*/
#include "NLRootMatrixSimpleFP.h"

void NLRootMatrixVerySimpleFPN10::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl) const {
    unsigned int num_atoms = 10;
    unsigned int max_distance = 3;
    addGenernalizedRepresentationFeature(fv, nl, num_atoms, max_distance);
}

void NLRootMatrixSimpleFPN8D3::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl) const {
    unsigned int num_atoms = 8;
    unsigned int max_distance = 3;
    bool include_adjacency_matrix = false;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, max_distance, include_adjacency_matrix);
}

void NLRootMatrixSimpleFPN10::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl) const {
    unsigned int num_atoms = 10;
    bool include_adjacency_matrix = false;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix);
}

void NLRootMatrixSimpleFPN16::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl) const {

    bool include_adjacency_matrix = false;
    unsigned int num_atoms = 16;
    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix);
}

void NLRootMatrixSimpleFPN32::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl) const {

    bool include_adjacency_matrix = false;
    unsigned int num_atoms = 32;
    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix);
}