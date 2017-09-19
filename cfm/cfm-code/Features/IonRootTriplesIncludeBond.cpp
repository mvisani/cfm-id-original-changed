/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootTriplesIncludeBond.cpp
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
#include "IonRootTriplesIncludeBond.h"

void IonRootTriplesIncludeBond::compute(FeatureVector &fv,
                                        const RootedROMolPtr *ion,
                                        const RootedROMolPtr *nl) const {

  int ring_break;
  const int path_len = 3;
  nl->mol.get()->getProp("IsRingBreak", ring_break);
  std::vector<path_t> paths;
  computeRootPaths(paths, ion, path_len, ring_break, true);
  addRootFeaturesWithBond(fv, paths, ring_break, path_len);
}