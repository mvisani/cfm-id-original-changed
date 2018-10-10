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

int GasteigerCharges::discretizeGasteigerCharge(double gc) const {

    //Discretize the Gasteiger Charges into 6 levels:
    // x < -0.5
    // -0.5 <= x < -0.1
    // -0.1 <= x < 0
    // 0 <= x < 0.1
    // 0.1 <= x <= 0.5
    // 0.5 <= x

    if (gc < -0.5) return 0;
    else if (gc < -0.1) return 1;
    else if (gc < 0) return 2;
    else if (gc < 0.1) return 3;
    else if (gc < 0.5) return 4;
    else return 5;
}

void
GasteigerCharges::compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl,
                          const int depth) const {

    //Collect the charges from the root atoms
    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);
    typedef std::pair<double, double> gasteiger_t;
    std::vector<gasteiger_t> gasteigers;

    //Ion
    double icharge, iothercharge;
    ion->root->getProp<double>("OrigGasteigerCharge", icharge);

    //Neutral Loss
    double nlcharge, nlothercharge;
    nl->root->getProp<double>("OrigGasteigerCharge", nlcharge);

    //Collate the charges
    gasteigers.push_back(gasteiger_t(icharge, nlcharge));
    if (ring_break) gasteigers.push_back(gasteiger_t(iothercharge, nlothercharge));

    if (!ring_break) {
        //Then there are 6 x 6 = 36 possible configurations of the charge
        // - Allocate one bit to each
        int gc_ion = discretizeGasteigerCharge(gasteigers[0].first);
        int gc_nl = discretizeGasteigerCharge(gasteigers[0].second);
        for (int i = 0; i <= 5; i++) {
            for (int j = 0; j <= 5; j++) {
                if (i == gc_ion && j == gc_nl) fv.addFeature(1.0);
                else fv.addFeature(0.0);
            }
        }
        //Ring break charges
        for (int i = 0; i < 36; i++)
            fv.addFeature(0.0);
    } else {
        //Non-Ring charges
        for (int i = 0; i < 36; i++)
            fv.addFeature(0.0);

        //Ring Charges - set bit if either breaks fit the rule
        int gc0_ion = discretizeGasteigerCharge(gasteigers[0].first);
        int gc0_nl = discretizeGasteigerCharge(gasteigers[0].second);
        int gc1_ion = discretizeGasteigerCharge(gasteigers[1].first);
        int gc1_nl = discretizeGasteigerCharge(gasteigers[1].second);

        for (int i = 0; i <= 5; i++) {
            for (int j = 0; j <= 5; j++) {
                if (i == gc0_ion && j == gc0_nl) fv.addFeature(1.0);
                else if (i == gc1_ion && j == gc1_nl) fv.addFeature(1.0);
                else fv.addFeature(0.0);
            }
        }
    }
}