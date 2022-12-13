/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentFingerPrintFeature.h
#
# Description: 	Classes for Fragment Finger Print Features
param.cpp.
#
# Copyright (c) 2013,2018
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FragmentFunctionalGroupFeature.h"

#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>


void FragmentFunctionalGroupFeature::compute(FeatureVector &fv, romol_ptr_t precursor_ion) const {

    auto fparams = new RDKit::FragCatParams(FGRPS_PICKLE);

    for(const auto & fg_param :  fparams->getFuncGroups()){

        std::vector<RDKit::MatchVectType> fgpMatches;
        int has_fg = RDKit::SubstructMatch(*precursor_ion, *fg_param,fgpMatches);
        if(has_fg){
            fv.addFeature(1.0);
        }
        else{
            fv.addFeature(0.0);
        }
    }

    delete fparams;
}