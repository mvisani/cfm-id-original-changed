/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# GraphDepthFeature.h
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

#ifndef CFM_GRAPHDEPTHFEATURE_H
#define CFM_GRAPHDEPTHFEATURE_H

#include "../Feature.h"

class GraphDepthFeature: public BreakFeature{
    public:
        GraphDepthFeature() {
            size = 5;
            name = "GraphDepth";
        };
        void
        compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const override;
};
#endif //CFM_GRAPHDEPTHFEATURE_H
