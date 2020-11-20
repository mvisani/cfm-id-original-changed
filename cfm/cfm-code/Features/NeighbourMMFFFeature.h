/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# NeighbourMMFFFeature.h
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
#pragma once

#include "../Feature.h"

class NeighbourMMFFFeature : public BreakFeature {
protected:
    void addNeighbourAtomTypes(FeatureVector &fv, const RootedROMol *mol,
                               const RDKit::Atom *root, int offset) const;
};