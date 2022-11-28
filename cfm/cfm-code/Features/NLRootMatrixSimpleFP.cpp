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

void NLRootGeneralizedMatrixFPN8::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {
    unsigned int num_atoms = 8;
    unsigned int max_distance = 8;
    addGenernalizedRepresentationFeature(fv, nl, num_atoms, max_distance);
}


void NLRootGeneralizedMatrixFPN10::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {
    unsigned int num_atoms = 10;
    unsigned int max_distance = 10;
    addGenernalizedRepresentationFeature(fv, nl, num_atoms, max_distance);
}

void NLRootMatrixSimpleFPN8D3::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {
    unsigned int num_atoms = 8;
    unsigned int max_distance = 3;
    bool include_adjacency_matrix = false;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, max_distance, include_adjacency_matrix, use_full_symbol_set);
}

void NLRootMatrixSimpleFPN10::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {
    unsigned int num_atoms = 10;
    bool include_adjacency_matrix = false;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}

void NLRootMatrixSimpleFPN16::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    bool include_adjacency_matrix = false;
    unsigned int num_atoms = 16;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}

void NLRootMatrixSimpleFPN32::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    bool include_adjacency_matrix = false;
    unsigned int num_atoms = 32;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}