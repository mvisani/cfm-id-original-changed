/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentFingerPrintFeature.cpp
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see
param.cpp.
#
# Copyright (c) 2013,2018
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include "FragmentFingerPrintFeature.h"

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/MolOps.h>

void FragmentFingerPrintFeature::compute(FeatureVector &fv, romol_ptr_t precursor_ion) const{

    int morgan_radius = 2;
    ExplicitBitVect *finger_print =
            RDKit::MorganFingerprints::getFingerprintAsBitVect(*precursor_ion, morgan_radius,
                                                               size);

    for (unsigned int i = 0; i < finger_print->getNumBits(); ++i)
        fv.addFeature((*finger_print)[i]);

    delete finger_print;
}