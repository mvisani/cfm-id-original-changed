
/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootEncodingD3Short.h
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
class IonRootEncodingD3Short : public FingerPirntFeature {
public:
  IonRootEncodingD3Short() {
    size = 512;
    name = "IonRootEncodingD3Short";
  };

  void compute(FeatureVector &fv, const RootedROMolPtr *ion,
               const RootedROMolPtr *nl) const;
};
>>>>>>> f8750844127ca8974966a0482f78619025ecfb54:cfm/cfm-code/Features/IonRootEncodingD3Short.h
