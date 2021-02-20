/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NeighbourMMFFFeature.cpp
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
#include "NeighbourMMFFFeature.h"

void NeighbourMMFFFeature::addNeighbourAtomTypes(FeatureVector &fv,
                                                 const RootedROMol *mol,
                                                 const RDKit::Atom *root,
                                                 int offset) const {
    // Iterate over the neighbours of the root atom
    RDKit::ROMol::ADJ_ITER_PAIR itp = mol->mol->getAtomNeighbors(root);
    int num_added = 0;
    for (; itp.first != itp.second; ++itp.first) {

        RDKit::Atom *nbr_atom = mol->mol->getAtomWithIdx(*itp.first);
        int atomtype;
        nbr_atom->getProp<int>("MMFFAtomType", atomtype);
        fv.addFeatureAtIdx(1.0, offset + atomtype);
        if (atomtype < 1 || atomtype > 99)
            fv.addFeatureAtIdx(1.0, offset + 100);
        num_added++;
    }
    // Additional feature indicating no neighbours
    if (num_added == 0)
        fv.addFeatureAtIdx(1.0, offset + 101);
}
