/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentFingerPrintFeature.h
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see
param.cpp.
#
# Copyright (c) 2013,2018
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#pragma once

#include "../Feature.h"

class FragmentFingerPrintFeature : public  FragmentFeature{
public:
    FragmentFingerPrintFeature() {
        size = 1024;
        name = "FragmentFingerPrint";
    };
    virtual void compute(FeatureVector &fv, romol_ptr_t precursor_ion) const override;
};

