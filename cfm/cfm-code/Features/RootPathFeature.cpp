/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# RootPathFeature.cpp
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
#include "RootPathFeature.h"

#include <GraphMol/MolOps.h>


void RootPathFeature::computeRootPaths(std::vector<path_t> &paths, const RootedROMol *mol, int len,
                                       bool with_bond = false) const {
    path_t path_so_far;
    addPathsFromAtom(paths, mol->root, mol->mol, mol->root, path_so_far, len,
                     with_bond);
}

void RootPathFeature::addPathsFromAtom(std::vector<path_t> &paths,
                                       const RDKit::Atom *atom,
                                       const romol_ptr_t mol,
                                       const RDKit::Atom *prev_atom,
                                       path_t &path_so_far, int len,
                                       bool with_bond) const {
    // Add the current symbol
    std::string symbol = atom->getSymbol();
    replaceUncommonWithX(symbol, false);
    path_so_far.push_back(symbol);

    // Iterate until len is reached, then add the path
    if (len > 1) {
        RDKit::ROMol::ADJ_ITER_PAIR itp = mol.get()->getAtomNeighbors(atom);
        for (; itp.first != itp.second; ++itp.first) {
            RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
            if (with_bond == true) {
                int bond_type = 0;
                mol.get()->getBondBetweenAtoms(atom->getIdx(), nbr_atom->getIdx())->getProp("OrigBondType", bond_type);
                path_so_far.push_back(boost::lexical_cast<std::string>(bond_type));
            }
            if (nbr_atom != prev_atom) {
                addPathsFromAtom(paths, nbr_atom, mol, atom, path_so_far, len - 1,
                                 with_bond);
            }
        }
    } else
        paths.push_back(path_so_far);

    // Remove the latest symbol
    path_so_far.pop_back();
}

void RootPathFeature::addRootPairFeatures(FeatureVector &fv,
                                          std::vector<path_t> &paths,
                                          int ring_break) const {
    // Add a feature indicating that there are no pairs
    fv.addFeature((double) (paths.size() == 0));

    // Iterate through all combinations of atom pairs, adding a count
    // of the number of times each is seen;
    // Note: the order matters here, root atom first then other
    std::vector<std::string>::const_iterator it1, it2;
    const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
    for (it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1) {
        for (it2 = ok_symbols->begin(); it2 != ok_symbols->end(); ++it2) {

            // Count how many of each possible symbol pair we have
            path_t sp;
            sp.push_back(*it1);
            sp.push_back(*it2);
            std::vector<path_t>::iterator it3 = paths.begin();
            double count = 0.0, ring_count = 0.0;
            for (; it3 != paths.end(); ++it3) {
                if (!ring_break && sp[0] == (*it3)[0] && sp[1] == (*it3)[1])
                    count += 1.0;
                if (ring_break && sp[0] == (*it3)[0] && sp[1] == (*it3)[1])
                    ring_count += 1.0;
            }

            // First feature indicates at least 1
            // Second feature indicates more than 1
            // Non-Ring
            if (count > 0.0)
                fv.addFeature(1.0);
            else
                fv.addFeature(0.0);
            if (count > 1.0)
                fv.addFeature(1.0);
            else
                fv.addFeature(0.0);
            // Ring
            if (ring_count > 0.0)
                fv.addFeature(1.0);
            else
                fv.addFeature(0.0);
            if (ring_count > 1.0)
                fv.addFeature(1.0);
            else
                fv.addFeature(0.0);
        }
    }
}

void RootPathFeature::addRootTripleFeatures(FeatureVector &fv,
                                            std::vector<path_t> &paths,
                                            int ring_break) const {
    // Add a feature indicating that there are no triples
    fv.addFeature((double) (paths.size() == 0));

    // Iterate through all combinations of atom triples, adding a count
    // of the number of times each is seen;
    // Note: the order matters here, root atom first then the next in the path,
    std::vector<std::string>::const_iterator it1, it2, it3;
    const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
    for (it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1) {
        for (it2 = ok_symbols->begin(); it2 != ok_symbols->end(); ++it2) {
            for (it3 = ok_symbols->begin(); it3 != ok_symbols->end(); ++it3) {
                // Count how many of each possible symbol pair we have
                path_t trp;
                trp.push_back(*it1);
                trp.push_back(*it2);
                trp.push_back(*it3);
                std::vector<path_t>::iterator it4 = paths.begin();
                double count = 0.0, ring_count = 0.0;
                for (; it4 != paths.end(); ++it4) {
                    if (trp[0] == (*it4)[0] && trp[1] == (*it4)[1] &&
                        trp[2] == (*it4)[2]) {
                        if (!ring_break)
                            count += 1.0;
                        else
                            ring_count += 1.0;
                    }
                }

                // First feature indicates at least 1
                // Second feature indicates more than 1
                // Non-Ring
                if (count > 0.0)
                    fv.addFeature(1.0);
                else
                    fv.addFeature(0.0);
                if (count > 1.0)
                    fv.addFeature(1.0);
                else
                    fv.addFeature(0.0);
                // Ring
                if (ring_count > 0.0)
                    fv.addFeature(1.0);
                else
                    fv.addFeature(0.0);
                if (ring_count > 1.0)
                    fv.addFeature(1.0);
                else
                    fv.addFeature(0.0);
            }
        }
    }
}

void RootPathFeature::addRootFeaturesWithBond(FeatureVector &fv,
                                              std::vector<path_t> &paths,
                                              int ring_break, int len) const {
    // double count = 0.0, ring_count = 0.0;
    int num_symbol_type = OKSymbolsLess().size();
    int base_feature_idx = fv.getTotalLength();
    std::unordered_map<int, int> counts;
    // for each path in the path list
    for (std::vector<path_t>::iterator path = paths.begin(); path != paths.end();
         ++path) {
        int feature_idx_offset = 0;
        for (unsigned int idx = 0; idx < path->size(); ++idx) {
            std::string symbol = path->at(idx);
            // if there are bond type information in the path
            // every second item in the list is bond type
            if (idx % 2 == 1) {
                feature_idx_offset =
                        num_symbol_type * feature_idx_offset + std::stoi(symbol);
            } else {
                feature_idx_offset =
                        num_symbol_type * feature_idx_offset + getSymbolsIndex(symbol, false);
            }
        }

        int key = ring_break ? feature_idx_offset * 4 : feature_idx_offset * 4 + 2;

        if (counts.find(key) == counts.end()) {
            counts[key] = 0;
        }
        counts[key] += 1;

        if (counts[key] > 1) {
            key += 1;
        }
        fv.addFeatureAtIdx(1.0, base_feature_idx + key);
    }
}
