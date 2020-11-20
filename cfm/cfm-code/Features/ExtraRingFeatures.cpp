/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# ExtraRingFeatures.cpp
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
#include "ExtraRingFeatures.h"
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>

void
ExtraRingFeatures::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    // Not a ring break
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);
    fv.addFeature(!ring_break);

    // Ion root is in ring
    // use findSSSR to init rinfo
    RDKit::MolOps::findSSSR(*ion->mol);
    RDKit::RingInfo *rinfo = ion->mol->getRingInfo();
    fv.addFeature(rinfo->minBondRingSize(ion->root->getIdx()) > 0);

    // NL root is in ring
    // use findSSSR to init rinfo
    RDKit::MolOps::findSSSR(*nl->mol);
    rinfo = nl->mol->getRingInfo();
    fv.addFeature(rinfo->minBondRingSize(nl->root->getIdx()) > 0);
}
