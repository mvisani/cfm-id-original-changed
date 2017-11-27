/*#########################################################################
# Mass Spec PredictNL and Identification of Metabolites
#
# NLRootEncodingD3.h
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

// Features use fingerprint encode NL fragmentatNL
class NLRootEncodingMorganD3 : public FingerPrintFeature {
public:
    NLRootEncodingMorganD3() {
        size = 1024;
        name = "NLRootEncodingMorganD3";
    };

    void compute(FeatureVector &fv, const RootedROMolPtr *NL,
                 const RootedROMolPtr *nl) const;
};
