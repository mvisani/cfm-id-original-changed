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
#include "FingerPrintFeature.h"
#pragma once

// Features use fingerprint encode NL
class IonRootGeneralizedMatrixFPN8 : public FingerPrintFeature {
public:
    IonRootGeneralizedMatrixFPN8() {
        size =  126 * 1 + 42 *7 ;//216 + 76 * 2; //48 + 28 ;
        name = "IonRootGeneralizedMatrixFPN8";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class IonRootGeneralizedMatrixFPN10 : public FingerPrintFeature {
public:
    IonRootGeneralizedMatrixFPN10() {
        size = 126 * 1 + 42 * 9; //60 + 50 + 45;
        name = "IonRootGeneralizedMatrixFPN10";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class IonRootMatrixSimpleFPN8D3 : public FingerPrintFeature {
public:
    IonRootMatrixSimpleFPN8D3() {
        size = 88;
        name = "IonRootMatrixSimpleFPN8D3";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class IonRootMatrixSimpleFPN10 : public FingerPrintFeature {
public:
    IonRootMatrixSimpleFPN10() {
        size = 110;
        name = "IonRootMatrixSimpleFPN10";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};

class IonRootMatrixSimpleFPN16 : public FingerPrintFeature {
public:
    IonRootMatrixSimpleFPN16() {
        size = 176;
        name = "IonRootMatrixSimpleFPN16";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};


class IonRootMatrixSimpleFPN32 : public FingerPrintFeature {
public:
    IonRootMatrixSimpleFPN32() {
        size = 352;
        name = "IonRootMatrixSimpleFPN32";
    };

    void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const override;
};