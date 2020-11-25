/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FunctionGroupFeature.h
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see param.cpp.
#
# Copyright (c) 2013,2017
# All rights reserved.
# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include "IonFunctionalGroupFeatures.h"

void IonFunctionalGroupFeatures::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);
    addFunctionalGroupFeatures(fv, ion, 1, false);
}

