/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# BreakHistoryFeature.cpp
#
# Description: 	Classes for BreakHistoryFeature
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include "BreakHistoryFeature.h"

void BreakHistoryFeature::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {


    std::string history_keyword = "FragmentationBondHistory";

    std::vector<int> history_so_far;
    ion->mol.get()->getProp(history_keyword, history_so_far);

    // we have seven bond types
    int occurs [7] = {0};

    // we only need the size -1 events, the latest event is current break
    // thus do not include it
    for(int i = 0 ; i < history_so_far.size() -1 ; ++i){
        // add occurs for each bond type
        occurs[history_so_far[i] - 1] += 1;
    }

    // This is a number
    // But for now let us trade it as categorical
    // TODO: visit this later
    for(auto & occur : occurs){
        for(int i = 1; i <= num_encoded_events; ++i){
            if(occur == i)
                fv.addFeature(1.0);
            else
                fv.addFeature(0.0);
        }
        if (occur > num_encoded_events)
            fv.addFeature(1.0);
        else
            fv.addFeature(0.0);
    }
}
