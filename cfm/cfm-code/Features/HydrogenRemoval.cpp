/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# HydrogenRemoval.h
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see param.cpp.
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include "HydrogenRemoval.h"

#include <GraphMol/PeriodicTable.h>
#include <GraphMol/AtomIterators.h>

void
HydrogenRemoval::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    double h_movement = 0.0;

    //Compute the mass difference in the neutral loss
    RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
    RDKit::ROMol::AtomIterator ai;
    for (ai = nl->mol.get()->beginAtoms(); ai != nl->mol.get()->endAtoms(); ++ai) {
        double orig_mass, mass = 0.0;
        std::string symbol = (*ai)->getSymbol();
        mass += pt->getMostCommonIsotopeMass(symbol);
        mass += (*ai)->getTotalNumHs() * pt->getMostCommonIsotopeMass("H");
        (*ai)->getProp<double>("OriginalMass", orig_mass);
        h_movement += (mass - orig_mass);
    }

    //Binary on/off indicating whether a particular transfer occurred
    for (double h = -4.0; h <= 4.0; h += 1.0) {
        if (fabs(h - h_movement) < 0.5)
            fv.addFeature(1.0);
        else
            fv.addFeature(0.0);
    }
    //Catch-all for all other hydrogen transfers
    if (fabs(h_movement) > 4.0)
        fv.addFeature(1.0);
    else
        fv.addFeature(0.0);
}
