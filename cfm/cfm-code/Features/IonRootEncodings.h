/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootEncodingD3.h
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

// Features use fingerprint encode ion fragmentation
class IonRootEncodingD3 : public FingerPrintFeature {
public:
    IonRootEncodingD3() {
        size = 512;
        name = "IonRootEncodingD3";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode ion fragmentation
class IonRootEncodingN10 : public FingerPrintFeature {
public:
    IonRootEncodingN10() {
        size = 1024;
        name = "IonRootEncodingN10";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode ion fragmentation
class IonRootEncodingD4 : public FingerPrintFeature {
public:
    IonRootEncodingD4() {
        size = 512;
        name = "IonRootEncodingD4";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode ion fragmentation
class IonRootEncodingMorganD3 : public FingerPrintFeature {
public:
    IonRootEncodingMorganD3() {
        size = 512;
        name = "IonRootEncodingMorganD3";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};