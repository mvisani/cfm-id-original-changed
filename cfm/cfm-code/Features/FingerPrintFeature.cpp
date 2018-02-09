/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FingerPrintFeature.h
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
#include "FingerPrintFeature.h"
#include "FeatureHelper.h"

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/MolOps.h>

#include <queue>

void FingerPrintFeature::getRemoveAtomIdx(
        const romol_ptr_t mol, const RDKit::Atom *root,
        std::vector<unsigned int> &remove_atom_ids, int distance) const {

    std::queue<const RDKit::Atom *> atom_queue;
    std::queue<int> distance_queue;
    atom_queue.push(root);
    distance_queue.push(0);

    std::unordered_set<unsigned int> visited;

    while (!atom_queue.empty() && !distance_queue.empty()) {
        const RDKit::Atom *curr = atom_queue.front();
        int curr_distance = distance_queue.front();
        atom_queue.pop();
        distance_queue.pop();

        // tracking cycles
        if (visited.find(curr->getIdx()) != visited.end()) {
            continue;
        }

        visited.insert(curr->getIdx());

        // thus we need remove it
        // also we does not care Hs
        if (distance < curr_distance || curr->getSymbol() == "H") {
            remove_atom_ids.push_back(curr->getIdx());
        }

        for (auto itp = mol->getAtomNeighbors(curr); itp.first != itp.second;
             ++itp.first) {
            const RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
            // if not parrent and I have not visit before
            if (curr != nbr_atom) {
                atom_queue.push(nbr_atom);
                distance_queue.push(curr_distance + 1);
            }
        }
    }
}

void FingerPrintFeature::removeAtomInTheList(
        RDKit::RWMol &mol, std::vector<unsigned int> &remove_atom_ids) const {

    // remove all duplications in the list
    // in theory this should not happen
    remove_atom_ids.erase(unique(remove_atom_ids.begin(), remove_atom_ids.end()),
                          remove_atom_ids.end());

    // Reverse sort list to make sure large idx get removed first
    // fail to do so will cause some problem since mol rerange idx after each
    // remove ops
    std::sort(remove_atom_ids.begin(), remove_atom_ids.end(),
              std::greater<int>());

    for (auto atom_idx : remove_atom_ids) {
        mol.removeAtom(atom_idx);
    }
}

// replace bond type with orig bond type
void FingerPrintFeature::replaceWithOrigBondType(RDKit::RWMol &rwmol) const {

    for (auto bi = rwmol.beginBonds(); bi != rwmol.endBonds(); ++bi) {
        int bond_type;
        (*bi)->getProp("OrigBondTypeRaw", bond_type);
        (*bi)->setBondType(static_cast<RDKit::Bond::BondType>(bond_type));
    }
}

void FingerPrintFeature::addRDKitFingerPrint(
        FeatureVector &fv, const RootedROMolPtr *mol, const RDKit::Atom *root,
        const unsigned int finger_print_size, const unsigned int max_nbr_distance,
        const unsigned int finger_print_min_path,
        const unsigned int finger_print_max_path) const {

    // Get list of atom we need to remove
    std::vector<unsigned int> remove_atom_ids;

    getRemoveAtomIdx(mol->mol, root, remove_atom_ids, max_nbr_distance);

    // Get Mol Object and remove atoms
    RDKit::RWMol part;
    part.insertMol(*(mol->mol));
    removeAtomInTheList(part, remove_atom_ids);

    // DO NOT sanitize mol because if we have a part of ring this going to break
    // RDKit::MolOps::sanitizeMol(part);

    // replace bond with OrigBondType
    // replaceWithOrigBondType(part);

    // Get finger prints with size
    ExplicitBitVect *fingerPrint = RDKit::RDKFingerprintMol(
            part, finger_print_min_path, finger_print_max_path, finger_print_size);

    for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
        fv.addFeature((*fingerPrint)[i]);
    }

    delete fingerPrint;
}

void FingerPrintFeature::addRDKitFingerPrintFeatures(
        FeatureVector &fv, const RootedROMolPtr *mol,
        unsigned int finger_print_size, unsigned int max_nbr_distance,
        int ring_break, unsigned int finger_print_min_path,
        unsigned int finger_print_max_path) const {

    addRDKitFingerPrint(fv, mol, mol->root, finger_print_size, max_nbr_distance,
                        finger_print_min_path, finger_print_max_path);

    if (ring_break) {
        addRDKitFingerPrint(fv, mol, mol->other_root, finger_print_size,
                            max_nbr_distance, finger_print_min_path,
                            finger_print_max_path);
    } else {
        for (int i = 0; i < finger_print_size; ++i) {
            fv.addFeature(0.0);
        }
    }
}

void FingerPrintFeature::addMorganFingerPrint(
        FeatureVector &fv, const RootedROMolPtr *mol, const RDKit::Atom *root,
        const unsigned int max_nbr_distance, const unsigned int finger_print_size,
        const int radius) const {

    // Get list of atom we need to remove
    std::vector<unsigned int> remove_atom_ids;
    std::unordered_set<unsigned int> visited;
    getRemoveAtomIdx(mol->mol, root, remove_atom_ids, max_nbr_distance);

    // Get Mol Object and remove atoms
    RDKit::RWMol part;
    part.insertMol(*(mol->mol));
    removeAtomInTheList(part, remove_atom_ids);

    // replace bond with OrigBondType
    // replaceWithOrigBondType(part);

    // Get finger prints with size
    ExplicitBitVect *fingerPrint =
            RDKit::MorganFingerprints::getFingerprintAsBitVect(part, radius,
                                                               finger_print_size);

    for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
        fv.addFeature((*fingerPrint)[i]);
    }

    delete fingerPrint;
}

void FingerPrintFeature::addMorganFingerPrintFeatures(
        FeatureVector &fv, const RootedROMolPtr *mol,
        const unsigned int finger_print_size, const unsigned int path_range,
        const int ring_break, const int radius) const {

    addMorganFingerPrint(fv, mol, mol->root, finger_print_size, path_range,
                         radius);

    if (ring_break) {
        addMorganFingerPrint(fv, mol, mol->other_root, finger_print_size,
                             path_range, radius);
    } else {
        for (int i = 0; i < finger_print_size; ++i) {
            fv.addFeature(0.0);
        }
    }
}

/*
std::string FingerPrintFeature::getSortingLabels(
    const romol_ptr_t mol, const RDKit::Atom *atom,
    const RDKit::Atom *prev_atom, int range,
    std::unordered_set<unsigned int> &visited,
    std::map<unsigned int, std::string> &sorting_labels) const {

  if (atom == nullptr || visited.find(atom->getIdx()) != visited.end()) {
    return "";
  }

  visited.insert(atom->getIdx());

  // get all the neighbors
  std::vector<std::string> children_keys;
  for (auto itp = mol->getAtomNeighbors(atom); itp.first != itp.second;
       ++itp.first) {
    RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
    if (nbr_atom != prev_atom) {
      std::string child_key;
      child_key = getSortingLabels(mol, nbr_atom, atom, range - 1, visited,
                                   sorting_labels);
      children_keys.push_back(child_key);
    }
  }

  // we don't care root node
  // there is no need to figure out which root to visit first anyway
  std::string atom_key = "";
  if (prev_atom != nullptr) {
    std::sort(children_keys.begin(), children_keys.end());
    // add bond type
    int bond_int = FeatureHelper::getBondTypeAsInt(
        mol->getBondBetweenAtoms(prev_atom->getIdx(), atom->getIdx()));
    atom_key += std::to_string(bond_int);

    // add atom symbol:
    // TODO: figure out which way is better
    // Maybe it is a better idea to use true symbol here
    // and replace uncommon with X later
    // replace symbol here or use true symbol
    std::string symbol_str = atom->getSymbol();
    // replaceUncommonWithX(symbol_str);
    atom_key += symbol_str;

    // get child atom keys str
    // sepeartor is important
    // otherwise CCC and C(C)C will has the same sorting key
    std::string children_atom_key = "|";
    for (auto child_key : children_keys) {
      children_atom_key += child_key;
    }
    sorting_labels.insert(std::pair<unsigned int, std::string>(
        atom->getIdx(), atom_key + children_atom_key));
  }
  return atom_key;
}
*/

std::string FingerPrintFeature::getSortingLabel(const romol_ptr_t mol,
                                                const RDKit::Atom *atom,
                                                const RDKit::Atom *parent_atom,
                                                int depth = 1) const {

    std::string atom_key;
    std::vector<std::string> children_keys;
    if (depth > 0) {
        for (auto itp = mol->getAtomNeighbors(atom); itp.first != itp.second;
             ++itp.first) {
            RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
            if (nbr_atom != parent_atom) {
                std::string child_key;
                child_key = getSortingLabel(mol, nbr_atom, atom, depth - 1);
                children_keys.push_back(child_key);
            }
        }
    }

    std::sort(children_keys.begin(), children_keys.end());
    // add bond type
    int bond_int = FeatureHelper::getBondTypeAsInt(
            mol->getBondBetweenAtoms(parent_atom->getIdx(), atom->getIdx()));
    atom_key += std::to_string(bond_int);

    std::string symbol_str = atom->getSymbol();
    // replaceUncommonWithX(symbol_str);
    if (symbol_str.size() == 1)
        symbol_str += " ";
    atom_key += symbol_str;

    // get child atom keys str
    std::string children_atom_key;
    for (const auto child_key : children_keys) {
        children_atom_key += child_key;
    }
    return atom_key;
}

// Method to get atom visited order via BFS
void FingerPrintFeature::getAtomVisitOrderBFS(
        const romol_ptr_t mol, const RDKit::Atom *root,
        std::vector<unsigned int> &visit_order, int range, int size,
        Bfs_Stop_Logic stop_logic) const {

    // maybe a struct is a better idea
    // but this is a one off
    std::queue<const RDKit::Atom *> atom_queue;
    std::queue<int> distance_queue;
    atom_queue.push(root);
    distance_queue.push(0);

    std::unordered_set<unsigned int> visited;

    while (!atom_queue.empty() && !distance_queue.empty()) {
        const RDKit::Atom *curr = atom_queue.front();
        int curr_distance = distance_queue.front();
        atom_queue.pop();
        distance_queue.pop();

        // if I have see this before
        if (std::find(visit_order.begin(), visit_order.end(), curr->getIdx()) !=
            visit_order.end()) {
            continue;
        }

        if (visit_order.size() < size) {
            visit_order.push_back(curr->getIdx());
        }

        if ((curr_distance < range && RANGE_ONLY == stop_logic) ||
            (visit_order.size() < size && SIZE_ONLY == stop_logic) ||
            (curr_distance < range && visit_order.size() < size &&
             RANGE_OR_SIZE == stop_logic)) {
            // use multimap since we can have duplicated labels
            std::multimap<std::string, const RDKit::Atom *> child_visit_order;

            for (auto itp = mol->getAtomNeighbors(curr); itp.first != itp.second;
                 ++itp.first) {
                RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
                // if we have not visit this node before
                // and this node is in the visit list
                if (nbr_atom != curr) {
                    std::string sorting_key = getSortingLabel(mol, nbr_atom, curr);
                    // std::string sorting_key = sorting_labels.at(nbr_atom->getIdx());
                    child_visit_order.insert(
                            std::pair<std::string, RDKit::Atom *>(sorting_key, nbr_atom));
                }
            }

            for (auto child : child_visit_order) {
                atom_queue.push(child.second);
                distance_queue.push(curr_distance + 1);
            }
        }
    }
}


void FingerPrintFeature::addAdjacentMatrixRepresentation(
        FeatureVector &fv, const RootedROMolPtr *mol, const RDKit::Atom *root,
        unsigned int max_nbr_distance, unsigned int num_atom, bool include_adjacency_matrix) const {

    // Get sorting labels
    // std::unordered_set<unsigned int> visited;
    // std::map<unsigned int, std::string> sorting_labels;
    // getSortingLabels(mol->mol, root, nullptr, num_atom, visited,
    // sorting_labels);

    // Get visit order
    std::vector<unsigned int> visit_order;
    getAtomVisitOrderBFS(mol->mol, root, visit_order, max_nbr_distance, num_atom,
                         SIZE_ONLY);

    std::map<unsigned int, int> visit_order_map;

    // make sure we only get num_atom amount of atoms, this is extra check, first check is done in the getAtomVisitOrderBFS
    for (int i = 0; i < visit_order.size() && i < num_atom; ++i) {
        visit_order_map[visit_order[i]] = i;
    }

    // init a 2D vector to store matrix
    std::vector<std::vector<int>> adjacency_matrix(num_atom,
                                                   std::vector<int>(num_atom, 0));

    // add bound type
    for (auto bi = mol->mol->beginBonds(); bi != mol->mol->endBonds(); ++bi) {
        // for each bond find two atoms
        unsigned int beginIdx = (*bi)->getBeginAtomIdx();
        unsigned int endIdx = (*bi)->getEndAtomIdx();
        int bond_type = FeatureHelper::getBondTypeAsInt(*bi);

        // if atoms in the list
        if (visit_order_map.find(beginIdx) != visit_order_map.end() &&
            visit_order_map.find(endIdx) != visit_order_map.end()) {
            adjacency_matrix[visit_order_map[beginIdx]][visit_order_map[endIdx]] =
                    bond_type;
            adjacency_matrix[visit_order_map[endIdx]][visit_order_map[beginIdx]] =
                    bond_type;
        }
    }

    if (true == include_adjacency_matrix) {
        // first bit indicate if there is a bond
        // rest 5 for each bond_type, one hot encoding
        const unsigned int num_bits_per_bond = 6;
        // we only need half of matrix exclude diagonal
        // so i start from 1 not 0
        for (int i = 0; i < num_atom; ++i) {
            for (int j = i + 1; j < num_atom; ++j) {
                int temp_feature[num_bits_per_bond] = {0};
                // check if bond exits/defined
                // std::cout << adjacency_matrix[i][j] << " ";
                if (adjacency_matrix[i][j] > 0) {
                    // first bit indicate if there is a bond
                    temp_feature[0] = 1;
                    // one hot encoding bond type
                    // since adjacency_matrix[i][j] > 0, bondtype 0 does not make sense
                    int bond_type = adjacency_matrix[i][j] > 5 ? 5 : adjacency_matrix[i][j];
                    temp_feature[bond_type] = 1;
                }
                // TODO Change to C++11 array
                fv.addFeatures(temp_feature, num_bits_per_bond);
            }
        }
    }

    // fv.printDebugInfo();
    // add atoms information into FP
    // TODO FIX THOSE MAGIC NUMBERS
    const unsigned int num_max_degree = 4;

    const unsigned int num_atom_types = 6;
    const unsigned int num_degree_feature_size = 5;
    for (int i = 0; i < num_atom; ++i) {
        int atom_type_feature[num_degree_feature_size] = {0};
        int atom_degree_feature[num_degree_feature_size] = {0};

        if (i < visit_order.size()) {
            int atom_idx = visit_order[i];

            // add atom types
            std::string symbol = mol->mol->getAtomWithIdx(atom_idx)->getSymbol();
            replaceUncommonWithX(symbol);
            int atom_feature = getSymbolsLessIndex(symbol);
            atom_type_feature[atom_feature] = 1;

            // add first order degree info
            int degree = mol->mol->getAtomWithIdx(atom_idx)->getDegree();
            degree = degree > num_max_degree ? num_max_degree : degree;
            atom_degree_feature[degree] = 1;

            // add second order degree info
            /*int second_order_;p[
            degree = mol->mol->getAtomWithIdx(atom_idx)->getDegree();
            degree = degree > num_max_degree ? num_max_degree : degree;
            atom_degree_feature[degree] = 1;*/
        }

        // TODO Change to C++11 array
        fv.addFeatures(atom_type_feature, num_atom_types);
        fv.addFeatures(atom_degree_feature, num_degree_feature_size);
    }
}

// for all the samples we have max atoms with a 3 atom group is 10
// for all the samples we have max atoms with a 5 atom group is 16
// therefore  we need 50 features for arcs
void FingerPrintFeature::addAdjacentMatrixRepresentationFeature(
        FeatureVector &fv, const RootedROMolPtr *mol,
        unsigned int max_nbr_distance, unsigned int num_atom,
        int ring_break, bool include_adjacency_matrix) const {

    addAdjacentMatrixRepresentation(fv, mol, mol->root, max_nbr_distance,
                                    num_atom, include_adjacency_matrix);
    if (ring_break > 0) {
        addAdjacentMatrixRepresentation(fv, mol, mol->other_root, max_nbr_distance,
                                        num_atom, include_adjacency_matrix);
    } else {
        // TODO: Get ride of this magic numbers
        unsigned int feature_size = num_atom * 11;
        if (true == include_adjacency_matrix) {
            feature_size += num_atom * (num_atom - 1) / 2 * 6;
        }
        for (int i = 0; i < feature_size; ++i) {
            fv.addFeature(0.0);
        }
    }
}