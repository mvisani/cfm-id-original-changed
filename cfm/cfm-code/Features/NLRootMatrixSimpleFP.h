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
#include "FingerPrintFeature.h"
#pragma once

// Features use fingerprint encode NL fragmentation
class NLRootGeneralizedMatrixFPN8 : public FingerPrintFeature {
public:
    NLRootGeneralizedMatrixFPN8() {
        size = 126 * 1 + 42 *7; // 60 + 45 + 50
        name = "NLRootGeneralizedMatrixFPN8";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootGeneralizedMatrixFPN10 : public FingerPrintFeature {
public:
    NLRootGeneralizedMatrixFPN10() {
        size =  126 * 1 + 42 * 9; // 60 + 45 + 50
        name = "NLRootGeneralizedMatrixFPN10";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixSimpleFPN8D3 : public FingerPrintFeature {
public:
    NLRootMatrixSimpleFPN8D3() {
        size = 88;
        name = "NLRootMatrixSimpleFPN8D3";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixSimpleFPN10 : public FingerPrintFeature {
public:
    NLRootMatrixSimpleFPN10() {
        size = 110;
        name = "NLRootMatrixSimpleFPN10";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixSimpleFPN16 : public FingerPrintFeature {
public:
    NLRootMatrixSimpleFPN16() {
        size = 176;
        name = "NLRootMatrixSimpleFPN16";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class NLRootMatrixSimpleFPN32 : public FingerPrintFeature {
public:
    NLRootMatrixSimpleFPN32() {
        size = 352;
        name = "NLRootMatrixSimpleFPN32";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};