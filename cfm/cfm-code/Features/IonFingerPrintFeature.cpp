/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonFingerPrintFeature.cpp
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
#include "IonFingerPrintFeature.h"

void
IonFingerPrintFeature::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    int morgan_radius = 2;
    addMorganFingerPrintFeatures(fv,ion,size, morgan_radius);
}