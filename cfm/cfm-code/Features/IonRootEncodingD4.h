/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootEncoding.h
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
#include "FingerPirntFeature.h"

// Features use fingerprint encode ion fragmentation
class IonRootEncodingD4 : public FingerPirntFeature {
public:
  IonRootEncodingD4() {
    size = 1024;
    name = "IonRootEncodingD4";
  };

  void compute(FeatureVector &fv, const RootedROMolPtr *ion,
               const RootedROMolPtr *nl) const;
};
