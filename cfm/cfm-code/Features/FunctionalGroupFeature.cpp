/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FunctionGroupFeature.cpp
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see param.cpp.
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include "FunctionalGroupFeature.h"

void FunctionalGroupFeature::addFunctionalGroupFeatures(FeatureVector &fv, const RootedROMol *mol, int max_depth,
                                                        bool extra) const {

    int offset = fv.getTotalLength();

    int num_grps = NUM_FGRPS;
    if (extra) num_grps = NUM_EXTRA_FGRPS;

    //Fill an array to begin with
    std::vector<int> tmp_full_fv((num_grps + 1) * (max_depth + 1), 0);
    addFunctionalGroupFeaturesFromAtom(tmp_full_fv, mol->root, mol->mol, mol->root, max_depth, 0, extra);

    //Then copy the array into the sparse format fv
    std::vector<int>::iterator it = tmp_full_fv.begin();
    for (int idx = 0; it != tmp_full_fv.end(); ++it, idx++)
        fv.addFeatureAtIdx((double) (*it), offset + idx);
}

void FunctionalGroupFeature::addFunctionalGroupFeaturesFromAtom(std::vector<int> &tmp_full_fv, const RDKit::Atom *atom,
                                                                const romol_ptr_t mol, const RDKit::Atom *prev_atom,
                                                                int max_depth, int depth, bool extra) const {

    int num_grps = NUM_FGRPS;
    if (extra)
        num_grps = NUM_EXTRA_FGRPS;

    //Check for functional groups at the current atom, and add them to the feature vector
    //iff they were not already found at a lesser depth.
    std::vector<unsigned int> fgrps;
    if (extra)
        atom->getProp<std::vector<unsigned int> >("ExtraFunctionalGroups", fgrps);
    else
        atom->getProp<std::vector<unsigned int> >("FunctionalGroups", fgrps);

    std::vector<unsigned int>::iterator it = fgrps.begin();
    for (; it != fgrps.end(); ++it) {
        bool added_at_lesser_depth = false;
        for (int d = 0; d <= depth; d++) {
            int idx = *it + d * (num_grps + 1);
            if (tmp_full_fv[idx]) {
                added_at_lesser_depth = true;
                break;
            }
        }
        if (!added_at_lesser_depth)
            tmp_full_fv[*it + depth * (num_grps + 1)] = 1;
    }

    //Iterate until max_depth is reached
    if (depth < max_depth) {
        RDKit::ROMol::ADJ_ITER_PAIR itp = mol.get()->getAtomNeighbors(atom);
        for (; itp.first != itp.second; ++itp.first) {
            RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
            if (nbr_atom != prev_atom)
                addFunctionalGroupFeaturesFromAtom(tmp_full_fv, nbr_atom, mol, atom, max_depth, depth + 1, extra);
        }
    }
}