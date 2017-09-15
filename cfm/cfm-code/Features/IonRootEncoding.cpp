/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootEncoding.cpp
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
#include "IonRootEncoding.h"
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/ExplicitBitVect.h>

void IonRootEncoding::compute(FeatureVector &fv,
                              const RootedROMolPtr *ion) const {
  RDKit::ROMol &ion_ref = *(ion->mol.get());

  unsigned int minPath = 1;
  unsigned int maxPath = 7;
  unsigned int fpSize = 2048;
  ExplicitBitVect *fingerPrint =
      RDKit::RDKFingerprintMol(ion_ref, minPath, maxPath, fpSize);
  for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
    fv.addFeature((*fingerPrint)[i]);
  }
}