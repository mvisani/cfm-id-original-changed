/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# BreakHistoryFeature.h
#
# Description: 	BreakHistoryFeature Class
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#pragma once

#include "../Feature.h"

class BreakHistoryFeature : public BreakFeature {
public:
    BreakHistoryFeature() {
        size = 7 * (num_encoded_events + 1); // 7 bond type and up to 6 for each
        name = "BreakHistoryFeature";
    };

    void
    compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const;

private:
    const int num_encoded_events = 5;
};