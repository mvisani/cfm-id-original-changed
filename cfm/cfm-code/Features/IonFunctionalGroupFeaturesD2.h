/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonFunctionalGroupFeaturesD2.h
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

#include "FunctionalGroupFeature.h"

class IonFunctionalGroupFeaturesD2 : public FunctionalGroupFeature {
public:
    IonFunctionalGroupFeaturesD2() {
        size = (NUM_FGRPS + 1) * 3;
        name = "IonFunctionalGroupFeaturesD2";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};
