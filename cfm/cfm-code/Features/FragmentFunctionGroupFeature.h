/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentFingerPrintFeature.h
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see
param.cpp.
#
# Copyright (c) 2013,2018
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/


#ifndef CFM_FRAGMENTFUNCTIONGROUP_H
#define CFM_FRAGMENTFUNCTIONGROUP_H

#include "../Feature.h"
#include "../FunctionalGroups.h"


class FragmentFunctionGroupFeature : public FragmentFeature {
public:
    FragmentFunctionGroupFeature() {
        size = NUM_FGRPS;
        name = "FragmentFunctionGroup";
    };
    virtual void compute(FeatureVector &fv, romol_ptr_t precursor_ion) const override;
};

#endif //CFM_FRAGMENTFUNCTIONGROUP_H
