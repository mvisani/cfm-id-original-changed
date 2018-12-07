/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# BrokenOrigBondType.cpp
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
#include "BrokenOrigBondType.h"

void BrokenOrigBondType::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl) const {
    int bondtype;
    nl->mol.get()->getProp("BrokenOrigBondType", bondtype);
    fv.addFeature(bondtype == 1); //SINGLE
    fv.addFeature(bondtype == 2); //DOUBLE
    fv.addFeature(bondtype == 3); //TRIPLE
    fv.addFeature(bondtype == 4); //AROMATIC
    fv.addFeature(bondtype == 5); //CONJUGATED
    fv.addFeature(bondtype == 6); //IONIC
    fv.addFeature(bondtype == 7); //H ONLY
}