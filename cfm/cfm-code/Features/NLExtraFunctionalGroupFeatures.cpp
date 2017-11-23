/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NLExtraFunctionalGroupFeatures.cpp
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
#include "NLExtraFunctionalGroupFeatures.h"

void NLExtraFunctionalGroupFeatures::compute(FeatureVector &fv,
                                             const RootedROMolPtr *ion,
                                             const RootedROMolPtr *nl) const {
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);
    addFunctionalGroupFeatures(fv, nl, 1, ring_break, true);
}