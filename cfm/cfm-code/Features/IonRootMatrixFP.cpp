/*#########################################################################
# Mass Spec PredictIon and Identification of Metabolites
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
##################################################### ####################*/
#include "IonRootMatrixFP.h"

void IonRootMatrixFPN6::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 6;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixFPN6D2::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 6;
    unsigned int max_distance = 2;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, max_distance, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixFPN8::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 8;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixFPN8D3::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 8;
    unsigned int max_distance = 3;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, max_distance, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixFPN10::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 10;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}

void IonRootMatrixFPN16::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 16;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}


void IonRootMatrixFPN10MoreSymbols::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 10;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = true;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}


void IonRootMatrixFPN16MoreSymbols::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int num_atoms = 16;
    bool include_adjacency_matrix = true;
    bool use_full_symbol_set = true;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, num_atoms, include_adjacency_matrix, use_full_symbol_set);
}