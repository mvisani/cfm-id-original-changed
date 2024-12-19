/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# GasteigerCharges.cpp
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
#include "GasteigerCharges.h"

int GasteigerCharges::discretizeGasteigerCharge(double gc) {

	// Discretize the Gasteiger Charges into 6 levels:
	//  x < -0.5
	//  -0.5 <= x < -0.1
	//  -0.1 <= x < 0
	//  0 <= x < 0.1
	//  0.1 <= x <= 0.5
	//  0.5 <= x

	if (gc < -0.5)
		return 0;
	else if (gc < -0.1)
		return 1;
	else if (gc < 0)
		return 2;
	else if (gc < 0.1)
		return 3;
	else if (gc < 0.5)
		return 4;
	else
		return 5;
}

void GasteigerCharges::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

	// Collect the charges from the root atoms
	// Ion
	double icharge;
	ion->root->getProp<double>("OrigGasteigerCharge", icharge);
	// Neutral Loss
	double nlcharge;
	nl->root->getProp<double>("OrigGasteigerCharge", nlcharge);

	// Collate the charges
	int gc_ion = discretizeGasteigerCharge(icharge);
	int gc_nl  = discretizeGasteigerCharge(nlcharge);
	for (int i = 0; i <= 5; i++) {
		for (int j = 0; j <= 5; j++) {
			if (i == gc_ion && j == gc_nl)
				fv.addFeature(1.0);
			else
				fv.addFeature(0.0);
		}
	}
}
