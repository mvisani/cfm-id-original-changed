/*#########################################################################
# Mass Spec PredictNL and Identification of Metabolites
#
# NLRootMatrixFP.h
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

// Features use fingerprint encode NL fragmentation
class NLRootMatrixFP : public FingerPrintFeature {
public:
    NLRootMatrixFP() {
        size = 380;
        name = "NLRootMatrixFP";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *NL,
                 const RootedROMolPtr *nl) const;
};
