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
#include "IonRootEncodings.h"

void IonRootEncodingD3::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int min_path = 1;
    unsigned int max_path = 2;
    unsigned int distance_to_root = 2;
    unsigned int finger_print_size = 512;

    addRDKitFingerPrintFeatures(fv, ion, finger_print_size, distance_to_root, true, min_path, max_path);
}

void IonRootEncodingN10::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int min_path = 1;
    unsigned int max_path = 3;
    unsigned int atom_count = 10;
    unsigned finger_print_size = 1024;

    addRDKitFingerPrintFeatures(fv, ion, finger_print_size, atom_count, true, min_path, max_path);
}

void IonRootEncodingD4::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int min_path = 1;
    unsigned int max_path = 3;
    unsigned int distance_to_root = 3;
    unsigned finger_print_size = 512;

    addRDKitFingerPrintFeatures(fv, ion, finger_print_size, distance_to_root, true, min_path, max_path);
}

void IonRootEncodingMorganD3::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    unsigned int distance_to_root = 2;
    unsigned int finger_print_size = 512;
    unsigned int morgan_radius = 2;

    addMorganFingerPrintFeatures(fv, ion, finger_print_size, distance_to_root, morgan_radius);
}