/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NLRootPairs.cpp
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
#include "NLRootPairs.h"

void NLRootPairs::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);
    std::vector<path_t> paths;
    computeRootPaths(paths, nl, 2, false);
    addRootPairFeatures(fv, paths, ring_break);
}
