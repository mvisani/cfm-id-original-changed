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

// feature_size = num_atom * 11  +   num_atom * (num_atom - 1) / 2 * 6;
// Features use fingerprint encode NL fragmentation

class NLRootMatrixFPN6 : public FingerPrintFeature {
public:
    NLRootMatrixFPN6() {
        size = 156;
        name = "NLRootMatrixFPN6";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixFPN6D2 : public FingerPrintFeature {
public:
    NLRootMatrixFPN6D2() {
        size = 156;
        name = "NLRootMatrixFPN6D2";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixFPN8 : public FingerPrintFeature {
public:
    NLRootMatrixFPN8() {
        size = 256;
        name = "NLRootMatrixFPN8";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixFPN8D3 : public FingerPrintFeature {
public:
    NLRootMatrixFPN8D3(){
        size = 256;
        name = "NLRootMatrixFPN8D3";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode NL fragmentation
class NLRootMatrixFPN10 : public FingerPrintFeature {
public:
    NLRootMatrixFPN10() {
        size = 380;
        name = "NLRootMatrixFPN10";
    };
    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixFPN16 : public FingerPrintFeature {
public:
    NLRootMatrixFPN16() {
        size = 896;
        name = "NLRootMatrixFPN16";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

// Features use fingerprint encode NL fragmentation with more symbols
class NLRootMatrixFPN10MoreSymbols : public FingerPrintFeature {
public:
    NLRootMatrixFPN10MoreSymbols() {
        size = 440 ;
        name = "NLRootMatrixFPN10MoreSymbols";
    };
    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixFPN16MoreSymbols : public FingerPrintFeature {
public:
    NLRootMatrixFPN16MoreSymbols() {
        size = 992;
        name = "NLRootMatrixFPN16MoreSymbols";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};