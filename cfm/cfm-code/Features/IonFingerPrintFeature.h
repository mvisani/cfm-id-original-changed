/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonFingerPrintFeature.h
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

#pragma once

#include "FingerPrintFeature.h"

class IonFingerPrintFeature : public FingerPrintFeature  {
public:
    IonFingerPrintFeature() {
        size =  512;
        name = "IonFingerPrintFeature";
    };

    void
    compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};
