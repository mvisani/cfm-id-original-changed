/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IonicFeatures.h
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
#include "IonicFeatures.h"
#include <GraphMol/AtomIterators.h>


void
IonicFeatures::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    int nl_pos = 0, nl_neg = 0, ion_pos = 0, ion_neg = 0;

    RDKit::ROMol::AtomIterator ai;
    for (ai = nl->mol.get()->beginAtoms(); ai != nl->mol.get()->endAtoms();
         ++ai) {
        int ionic_frag_q;
        (*ai)->getProp("IonicFragmentCharge", ionic_frag_q);
        if (ionic_frag_q < 0)
            nl_neg = 1;
        if (ionic_frag_q > 0)
            nl_pos = 1;
    }
    for (ai = ion->mol.get()->beginAtoms(); ai != ion->mol.get()->endAtoms();
         ++ai) {
        int ionic_frag_q;
        (*ai)->getProp("IonicFragmentCharge", ionic_frag_q);
        if (ionic_frag_q < 0)
            ion_neg = 1;
        if (ionic_frag_q > 0)
            ion_pos = 1;
    }

    fv.addFeature(nl_pos);  // NL has positive ionic fragment
    fv.addFeature(ion_pos); // Ion has positive ionic fragment
    fv.addFeature(nl_neg);  // NL has negative ionic fragment
    fv.addFeature(ion_neg); // Ion has negative ionic fragment
    fv.addFeature(!(nl_pos || nl_neg || ion_neg || ion_pos)); // No ionic fragments anywhere
}
