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
#########################################################################*/
#include "NLRootMatrixSimpleFP.h"

void NLRootMatrixSimpleFP::compute(FeatureVector &fv, const RootedROMolPtr *ion,
                                   const RootedROMolPtr *nl, const int depth) const {
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);

    unsigned int num_atoms = 10;
    bool include_adjacency_matrix = false;

    addAdjacentMatrixRepresentationFeature(fv, nl, num_atoms, ring_break, include_adjacency_matrix);
}