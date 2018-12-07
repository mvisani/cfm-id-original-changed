/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# GraphDepthFeature.cpp
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
#include "GraphDepthFeature.h"

void GraphDepthFeature::compute(FeatureVector &fv, romol_ptr_t precursor_ion, int depth) const {
    std::vector<int> features(this->size,0);
    if(depth < features.size())
        features[depth] = 1;
    fv.addFeatures(features);
}