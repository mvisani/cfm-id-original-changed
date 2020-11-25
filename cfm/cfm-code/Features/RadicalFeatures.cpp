/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# RadicalFeatures.cpp
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
#include "RadicalFeatures.h"

void
RadicalFeatures::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    int ion_radical = moleculeHasSingleRadical(ion->mol.get());
    int nl_radical = moleculeHasSingleRadical(nl->mol.get());
    fv.addFeature(ion_radical);                 // Ion is radical
    fv.addFeature(nl_radical);                  // NL is radical
    fv.addFeature(!ion_radical && !nl_radical); // Neither NL or Ion are radical
}
