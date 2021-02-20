/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FunctionGroupFeature.h
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
#include "../FunctionalGroups.h"

class FunctionalGroupFeature : public BreakFeature {
protected:
    void addFunctionalGroupFeaturesFromAtom(std::vector<int> &tmp_full_fv,
                                            const RDKit::Atom *atom,
                                            romol_ptr_t mol,
                                            const RDKit::Atom *prev_atom,
                                            int max_depth, int depth,
                                            bool extra) const;

    void addFunctionalGroupFeatures(FeatureVector &fv, const RootedROMol *mol, int max_dist,
                                        bool extra) const;
};
