/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonRootMMFFAtomType.cpp
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
#include "IonRootMMFFAtomType.h"

void IonRootMMFFAtomType::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {
    int offset = fv.getTotalLength() - 1;
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);

    // Set features for atom types
    int atomtype, otheratomtype = 0;
    ion->root->getProp<int>("MMFFAtomType", atomtype);
    fv.addFeatureAtIdx(1.0, offset + atomtype);

    // 100 Features in total - last features indicates out-of-range
    if (atomtype < 1 || atomtype > 99)
        fv.addFeatureAtIdx(1.0, offset + 100);
    else if (ring_break && (otheratomtype < 1 || otheratomtype > 99))
        fv.addFeatureAtIdx(1.0, offset + 100);
    else
        fv.addFeatureAtIdx(0.0, offset + 100);
}