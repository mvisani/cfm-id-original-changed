/*#########################################################################
# Mass Spec PredictNL and Identification of Metabolites
#
# IonRootMatrixFP.cpp
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
#########################################################################*/
#include "IonRootMatrixSimpleFP.h"

void IonRootGeneralizedMatrixFPN8::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 8;
    unsigned int max_distance = 8;
    addGenernalizedRepresentationFeature(fv, ion, num_atoms, max_distance);
}

void IonRootGeneralizedMatrixFPN10::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 10;
    unsigned int max_distance = 10;
    addGenernalizedRepresentationFeature(fv, ion, num_atoms, max_distance);
}

void IonRootMatrixSimpleFPN8D3::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 8;
    unsigned int max_distance = 3;
    bool include_adjacency_matrix = false;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, max_distance, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixSimpleFPN10::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 10;
    bool include_adjacency_matrix = false;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixSimpleFPN16::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 16;
    bool include_adjacency_matrix = false;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixSimpleFPN32::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 32;
    bool include_adjacency_matrix = false;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}