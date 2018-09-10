/*#########################################################################
# Mass Spec PredictNL and Identification of Metabolites
#
# IonRootMatrixFP.h
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
#pragma once
#include "FingerPrintFeature.h"

// Features use fingerprint encode NL fragmentatNL

//feature_size = num_atom * 11 + 5 * (num_atom) +   num_atom * (num_atom - 1) / 2 * 6;
// Features use fingerprint encode NL fragmentation
class IonRootMatrixFPN6 : public FingerPrintFeature {
public:
    IonRootMatrixFPN6() {
        size = 372;
        name = "IonRootMatrixFPN6";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};

class IonRootMatrixFPN8 : public FingerPrintFeature {
public:
    IonRootMatrixFPN8() {
        size = 592;
        name = "IonRootMatrixFPN8";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};

class IonRootMatrixFPN10 : public FingerPrintFeature {
public:
    IonRootMatrixFPN10() {
        size = 860;
        name = "IonRootMatrixFPN10";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};

class IonRootMatrixFPN16 : public FingerPrintFeature {
public:
    IonRootMatrixFPN16() {
        size = 1952;
        name = "IonRootMatrixFPN16";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};