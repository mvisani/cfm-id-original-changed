/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootEncodingD3.cpp
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
#include "IonRootEncodingMorganD3.h"

void IonRootEncodingMorganD3::compute(FeatureVector &fv, const RootedROMolPtr *ion,
                                      const RootedROMolPtr *nl) const {
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);

    unsigned int path_range = 3;
    unsigned int finger_print_size = 512;
    unsigned int morgan_radius = 2;

    addMorganFingerPrint(fv, ion, finger_print_size, path_range, ring_break, morgan_radius);
}