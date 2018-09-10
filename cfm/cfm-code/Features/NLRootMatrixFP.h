/*#########################################################################
# Mass Spec PredictNL and Identification of Metabolites
#
# NLRootMatrixFPN10.h
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

//feature_size = num_atom * 11 + 5 * (num_atom) +   num_atom * (num_atom - 1) / 2 * 6;
// Features use fingerprint encode NL fragmentation

class NLRootMatrixFPN6 : public FingerPrintFeature {
public:
    NLRootMatrixFPN6() {
        size = 372;
        name = "NLRootMatrixFPN6";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};

class NLRootMatrixFPN8 : public FingerPrintFeature {
public:
    NLRootMatrixFPN8() {
        size = 592;
        name = "NLRootMatrixFPN8";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};

// Features use fingerprint encode NL fragmentation
class NLRootMatrixFPN10 : public FingerPrintFeature {
public:
    NLRootMatrixFPN10() {
        size = 860;
        name = "NLRootMatrixFPN10";
    };
    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};

class NLRootMatrixFPN16 : public FingerPrintFeature {
public:
    NLRootMatrixFPN16() {
        size = 1952;
        name = "NLRootMatrixFPN16";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};