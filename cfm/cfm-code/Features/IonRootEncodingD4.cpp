/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootEncoding.cpp
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
#include "IonRootEncodingD4.h"

void IonRootEncodingD4::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl,
                                const int depth) const {
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);

    unsigned int min_path = 1;
    unsigned int max_path = 3;
    unsigned int distance_to_root = 3;
    unsigned finger_print_size = 512;

    addRDKitFingerPrintFeatures(fv, ion, finger_print_size, distance_to_root, ring_break, true, min_path, max_path);
}