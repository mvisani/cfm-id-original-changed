/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NLRootEncodingD3.cpp
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
class NLRootEncodingD3 : public FingerPrintFeature {

public:
    NLRootEncodingD3() {
        size = 512;
        name = "NLRootEncodingD3";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode ion fragmentation
class NLRootEncodingN10 : public FingerPrintFeature {
public:
    NLRootEncodingN10() {
        size = 1024;
        name = "NLRootEncodingN10";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode NL fragmentatNL
class NLRootEncodingMorganD3 : public FingerPrintFeature {
public:
    NLRootEncodingMorganD3() {
        size = 512;
        name = "NLRootEncodingMorganD3";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode ion fragmentation
class NLRootEncodingD4 : public FingerPrintFeature {
public:
    NLRootEncodingD4() {
        size = 512;
        name = "NLRootEncodingD4";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};