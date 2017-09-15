/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FunctionGroupFeature.h
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
#include "IonNeighbourMMFFAtomType.h"

void IonNeighbourMMFFAtomType::compute(FeatureVector &fv,
                                       const RootedROMolPtr *ion,
                                       const RootedROMolPtr *nl) const {
  int offset = fv.getTotalLength() - 1;
  int ring_break;
  nl->mol.get()->getProp("IsRingBreak", ring_break);
  fv.addFeatureAtIdx(0.0,
                     offset + 101); // Make the feature vector the right length
  addNeighbourAtomTypes(fv, ion, ion->root, offset);
  if (ring_break)
    addNeighbourAtomTypes(fv, ion, ion->other_root, offset);
}
