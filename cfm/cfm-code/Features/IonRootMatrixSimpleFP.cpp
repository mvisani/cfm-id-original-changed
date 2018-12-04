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


void IonRootMatrixVerySimpleFPN10::compute(FeatureVector &fv, const RootedROMolPtr *ion,
                                       const RootedROMolPtr *nl, int depth) const {

    unsigned int num_atoms = 10;
    addGenernalizedRepresentationFeature(fv, ion, num_atoms);
}


void IonRootMatrixSimpleFPN10::compute(FeatureVector &fv, const RootedROMolPtr *ion,
                                    const RootedROMolPtr *nl, int depth) const {

    unsigned int num_atoms = 10;
    bool include_adjacency_matrix = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, 0, include_adjacency_matrix);
}

void IonRootMatrixSimpleFPN16::compute(FeatureVector &fv, const RootedROMolPtr *ion,
                                    const RootedROMolPtr *nl, int depth) const {

    unsigned int num_atoms = 16;
    bool include_adjacency_matrix = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, 0, include_adjacency_matrix);
}

void IonRootMatrixSimpleFPN32::compute(FeatureVector &fv, const RootedROMolPtr *ion,
                                       const RootedROMolPtr *nl, int depth) const {

    unsigned int num_atoms = 32;
    bool include_adjacency_matrix = false;

    addAdjacentMatrixRepresentationFeature(fv, ion, num_atoms, 0, include_adjacency_matrix);
}