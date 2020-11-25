/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# RootAtomFeature.cpp
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
#include "RootAtomFeature.h"

void RootAtomFeature::computeRootAtomFeature(FeatureVector &fv,
                                             const RootedROMol *mol,
                                             bool ring_break) const {

    int offset = fv.getTotalLength();

    // Add the root atom
    std::vector<std::string>::const_iterator it1;
    const std::vector<std::string> *ok_symbols = &OKsymbols();
    int idx = 0;
    bool found = false;
    for (it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1, idx++) {
        if ((*it1) == mol->root->getSymbol()) {
            fv.addFeatureAtIdx(1.0, offset + idx);
            found = true;
        } else
            fv.addFeatureAtIdx(0.0, offset + idx);
    }
    if (mol->root->getSymbol() == "H") {
        fv.addFeatureAtIdx(1.0, offset + idx);
        found = true;
    } else
        fv.addFeatureAtIdx(0.0, offset + idx);
    idx++;
    if (!found)
        fv.addFeatureAtIdx(1.0, offset + idx);
    else
        fv.addFeatureAtIdx(0.0, offset + idx);
}