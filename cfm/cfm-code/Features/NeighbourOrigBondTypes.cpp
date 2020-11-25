/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NeighbourOrigBondTypes.cpp
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
#include "NeighbourOrigBondTypes.h"

void NeighbourOrigBondTypes::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);
    addNeighbourOrigBondFeatures(fv, ion, ring_break);
    addNeighbourOrigBondFeatures(fv, nl, ring_break);
}

void addNeighbourOrigBondFeatures(FeatureVector &fv, const RootedROMol *mol,
                                  int ring_break) {

    std::vector<int> seen_types(6, 0);
    int feature_offset = fv.getTotalLength();
    //RDKit::ROMol::ADJ_ITER_PAIR
    auto itp = mol->mol->getAtomNeighbors(mol->root);
    for (; itp.first != itp.second; ++itp.first) {
        RDKit::Bond *bond =
                mol->mol->getBondBetweenAtoms(*itp.first, mol->root->getIdx());
        int bondtype;
        bond->getProp("OrigBondType", bondtype);
        int idx = feature_offset + bondtype;
        if (!seen_types[bondtype])
            fv.addFeatureAtIdx(1.0, idx);
        seen_types[bondtype] = 1;
    }

    if (fv.getTotalLength() - feature_offset == 0)
        fv.addFeatureAtIdx(1.0, feature_offset); // No connected bonds
    if (fv.getTotalLength() - feature_offset < 6)
        fv.addFeatureAtIdx(0.0, feature_offset + 5); // Update length
}
