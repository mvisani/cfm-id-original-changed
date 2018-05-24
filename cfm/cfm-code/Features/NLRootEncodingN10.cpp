/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NLRootEncoding.cpp
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
#include "NLRootEncodingN10.h"

void NLRootEncodingN10::compute(FeatureVector &fv, const RootedROMolPtr *NL,
                                   const RootedROMolPtr *nl) const {
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);

    unsigned int min_path = 1;
    unsigned int max_path = 3;
    unsigned int atom_count = 10;
    unsigned finger_print_size = 512;

    addRDKitFingerPrintFeatures(fv, NL, finger_print_size, atom_count, ring_break, false, min_path, max_path);
}