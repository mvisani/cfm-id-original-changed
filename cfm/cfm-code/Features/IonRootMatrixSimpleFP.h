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
class IonRootMatrixSimpleFPN10 : public FingerPrintFeature {
public:
    IonRootMatrixSimpleFPN10() {
        size = 110;
        name = "IonRootMatrixSimpleFPN10";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};

class IonRootMatrixSimpleFPN16 : public FingerPrintFeature {
public:
    IonRootMatrixSimpleFPN16() {
        size = 176;
        name = "IonRootMatrixSimpleFPN16";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};


class IonRootMatrixSimpleFPN32 : public FingerPrintFeature {
public:
    IonRootMatrixSimpleFPN32() {
        size = 352;
        name = "IonRootMatrixSimpleFPN32";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};