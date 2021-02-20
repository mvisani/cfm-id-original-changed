/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# QuadraticFeatures.cpp
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

#include "QuadraticFeatures.h"

void
QuadraticFeatures::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {
    // Compute quadratic feature indexes for all existing features
    int n = fv.getTotalLength();
    std::vector<int> quadratic_indexes;
    for (auto it1 = std::next(fv.getFeatureBegin()); it1 != fv.getFeatureEnd(); ++it1) {
        // Due to symmetry and not wanting to include square features,
        // or bias features, only the lower left triangle of the
        // feature x feature matrix is included (minus bias row/col)
        //- the offset gives the index for the first used feature in each row.
        int offset = n + (*it1 - 2) * (*it1 - 1) / 2;
        for (auto it2 = std::next(fv.getFeatureBegin()); it2 != it1; ++it2)
            quadratic_indexes.push_back(offset + *it2 - 1);
    }
    // Add the features
    // Note: modifying the feature vector in the above loop causes problems...
    std::vector<int>::iterator it = quadratic_indexes.begin();
    for (; it != quadratic_indexes.end(); ++it)
        fv.addFeatureAtIdx(1.0, *it);

    // Update the feature length (without modifying any features).
    int total_num_features = n + (n - 1) * (n - 2) / 2;
    fv.addFeatureAtIdx(0.0, total_num_features - 1);
}