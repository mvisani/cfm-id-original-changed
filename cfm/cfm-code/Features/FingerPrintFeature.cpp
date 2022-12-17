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
#include <GraphMol/AtomIterators.h>

#include <queue>
#include <bitset>

void FingerPrintFeature::getRemoveAtomIdxByDisatnce(
        romol_ptr_t mol, const RDKit::Atom *root,
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

void FingerPrintFeature::getRemoveAtomIdxByCount(romol_ptr_t mol, const RDKit::Atom *root,
                                                 std::vector<unsigned int> &remove_atom_ids,
                                                 int count) const {
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
        if (visited.size() >= count) {
            remove_atom_ids.push_back(curr->getIdx());
        }

        // use multimap since we can have duplicated labels
        std::multimap<std::string, const RDKit::Atom *> child_visit_order;

        for (auto itp = mol->getAtomNeighbors(curr); itp.first != itp.second; ++itp.first) {
            RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
            // if we have not visit this node before
            // and this node is in the visit list
            if (nbr_atom != curr) {
                /*std::string sorting_key = getSortingLabel(mol, nbr_atom, curr, <#initializer#>);
                child_visit_order.insert(
                        std::pair<std::string, RDKit::Atom *>(sorting_key, nbr_atom));*/
            }
        }

        for (auto child : child_visit_order) {
            atom_queue.push(child.second);
            distance_queue.push(curr_distance + 1);
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

void FingerPrintFeature::addRDKitFingerPrint(std::vector<int> &tmp_fv, const RootedROMol *mol,
                                             const RDKit::Atom *root,
                                             unsigned int finger_print_size, unsigned int limitation_param,
                                             unsigned int finger_print_min_path, unsigned int finger_print_max_path,
                                             bool limited_by_distance) const {

    // Get list of atom we need to remove
    std::vector<unsigned int> remove_atom_ids;
    if (limited_by_distance)
        getRemoveAtomIdxByDisatnce(mol->mol, root, remove_atom_ids, limitation_param);
    else
        getRemoveAtomIdxByCount(mol->mol, root, remove_atom_ids, limitation_param);

    // Get Mol Object and remove atoms
    RDKit::RWMol part;
    part.insertMol(*(mol->mol));
    removeAtomInTheList(part, remove_atom_ids);

    // Get finger prints with size
    ExplicitBitVect *finger_print = RDKit::RDKFingerprintMol(
            part, finger_print_min_path, finger_print_max_path, finger_print_size);

    tmp_fv.resize(finger_print->getNumBits());
    for (unsigned int i = 0; i < finger_print->getNumBits(); ++i)
        tmp_fv[i] = (*finger_print)[i];

    delete finger_print;
}

void FingerPrintFeature::addRDKitFingerPrintFeatures(FeatureVector &fv, const RootedROMol *mol,
                                                     unsigned int finger_print_size,
                                                     unsigned int limitation_param, bool limited_by_distance,
                                                     unsigned int finger_print_min_path,
                                                     unsigned int finger_print_max_path) const {

    std::vector<int> local_tmp_fv;
    addRDKitFingerPrint(local_tmp_fv, mol, mol->root, finger_print_size, limitation_param, finger_print_min_path,
                        finger_print_max_path, limited_by_distance);
    fv.addFeatures(local_tmp_fv);

}

void FingerPrintFeature::addMorganFingerPrint(
        std::vector<int> &tmp_fv, const RootedROMol *mol, const RDKit::Atom *root,
        const unsigned int max_nbr_distance, const unsigned int finger_print_size,
        const int radius) const {

    // Get list of atom we need to remove
    std::vector<unsigned int> remove_atom_ids;
    std::unordered_set<unsigned int> visited;
    getRemoveAtomIdxByDisatnce(mol->mol, root, remove_atom_ids, max_nbr_distance);

    // Get Mol Object and remove atoms
    RDKit::RWMol part;
    part.insertMol(*(mol->mol));
    removeAtomInTheList(part, remove_atom_ids);

    // Get finger prints with size
    ExplicitBitVect *finger_print =
            RDKit::MorganFingerprints::getFingerprintAsBitVect(part, radius,
                                                               finger_print_size);

    tmp_fv.resize(finger_print->getNumBits());
    for (unsigned int i = 0; i < finger_print->getNumBits(); ++i)
        tmp_fv[i] = (*finger_print)[i];

    delete finger_print;
}

void FingerPrintFeature::addMorganFingerPrintFeatures(FeatureVector &fv, const RootedROMol *mol,
                                                      unsigned int finger_print_size,
                                                      unsigned int path_range, int radius) const {

    std::vector<int> local_tmp_fv;
    addMorganFingerPrint(local_tmp_fv, mol, mol->root, finger_print_size, path_range,
                         radius);

    fv.addFeatures(local_tmp_fv);
}


void FingerPrintFeature::addMorganFingerPrintFeatures(FeatureVector &fv,
                                                      const RootedROMol *mol, unsigned int finger_print_size,
                                                      int radius) const {

    std::vector<int> local_tmp_fv;

    // Get finger prints with size
    ExplicitBitVect *finger_print =
            RDKit::MorganFingerprints::getFingerprintAsBitVect((*mol->mol), radius,
                                                               finger_print_size);

    local_tmp_fv.resize(finger_print->getNumBits());
    for (unsigned int i = 0; i < finger_print->getNumBits(); ++i)
        local_tmp_fv[i] = (*finger_print)[i];


    fv.addFeatures(local_tmp_fv);
    delete finger_print;
}

std::string FingerPrintFeature::getSortingLabel(const romol_ptr_t mol, const RDKit::Atom *atom,
                                                std::map<unsigned int, unsigned int> &distances,
                                                std::map<unsigned int, std::string> &labels,
                                                bool use_full_symbols_set) const {
    auto distance_to_root = distances[atom->getIdx()];

    // add bond and atom label
    std::string atom_key = "";

    // in case we have more than one bond lead to this node
    if (distance_to_root > 0) {
        int bond_type = 1;
        for (auto itp = mol->getAtomNeighbors(atom); itp.first != itp.second; ++itp.first) {
            RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
            auto child_distance_to_root = distances[nbr_atom->getIdx()];
            int bond_int = FeatureHelper::getBondTypeAsInt(
                    mol->getBondBetweenAtoms(atom->getIdx(), nbr_atom->getIdx()));

            // use smallest bond order
            if (child_distance_to_root < distance_to_root && bond_int < bond_type)
                bond_type = bond_int;
        }
        atom_key += std::to_string(bond_type);
    }

    std::string symbol_str = atom->getSymbol();
    replaceUncommonWithX(symbol_str, use_full_symbols_set);
    if (symbol_str.size() == 1)
        symbol_str += " ";
    atom_key += symbol_str;

    std::vector<std::string> children_keys;
    for (auto itp = mol->getAtomNeighbors(atom); itp.first != itp.second; ++itp.first) {
        RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
        auto child_distance_to_root = distances[nbr_atom->getIdx()];
        // make sure we are not going back and not going in a loop
        if (child_distance_to_root > distance_to_root) {
            std::string child_key;
            child_key = getSortingLabel(mol, nbr_atom, distances, labels, use_full_symbols_set);
            children_keys.push_back(child_key);
        }
    }

    std::sort(children_keys.begin(), children_keys.end());

    // get child atom keys str
    std::string children_atom_key;
    for (const auto child_key : children_keys)
        children_atom_key +=  "[" + child_key +"]";
    // save the label
    if(!children_keys.empty())
        labels[atom->getIdx()] = atom_key + children_atom_key;
    else
        labels[atom->getIdx()] = atom_key;

    return atom_key;
}

void FingerPrintFeature::getAtomDistanceToRoot(const RootedROMol *roMolPtr,
                                                std::map<unsigned int, unsigned int> &distances) const {
    // maybe a struct is a better idea
    // but this is a one off
    std::queue<const RDKit::Atom *> atom_queue;
    std::queue<unsigned int> distance_queue;
    atom_queue.push(roMolPtr->root);
    distance_queue.push(0);

    auto mol = roMolPtr->mol;

    while (!atom_queue.empty()) {
        const RDKit::Atom *curr = atom_queue.front();
        auto curr_distance = distance_queue.front();

        atom_queue.pop();
        distance_queue.pop();

        // if I have see this before
        if (distances.find(curr->getIdx()) != distances.end())
            continue;
        distances[curr->getIdx()] = curr_distance;

        for (auto itp = mol->getAtomNeighbors(curr); itp.first != itp.second; ++itp.first) {
            RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
            atom_queue.push(nbr_atom);
            distance_queue.push(curr_distance + 1);
        }
    }
}

void FingerPrintFeature::getBondAtomPairAtEachDistance(const RootedROMol *roMolPtr,
                                                       std::vector<std::map<std::string, int>> &dict,
                                                       bool use_full_symbols_set) const {
    // maybe a struct is a better idea
    // but this is a one off
    std::queue<const RDKit::Atom *> atom_queue;
    std::queue<unsigned int> distance_queue;
    std::map<unsigned int, unsigned int> distances;
    atom_queue.push(roMolPtr->root);
    distance_queue.push(0);

    auto mol = roMolPtr->mol;

    while (!atom_queue.empty()) {
        const RDKit::Atom *curr = atom_queue.front();
        auto curr_distance = distance_queue.front();

        atom_queue.pop();
        distance_queue.pop();

        // if I have see this before
        if (distances.find(curr->getIdx()) != distances.end())
            continue;
        distances[curr->getIdx()] = curr_distance;

        for (auto itp = mol->getAtomNeighbors(curr); itp.first != itp.second; ++itp.first) {
            RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
            atom_queue.push(nbr_atom);
            distance_queue.push(curr_distance + 1);

            auto bi = mol.get()->getBondBetweenAtoms(curr->getIdx(), nbr_atom->getIdx());
            int bond_type = FeatureHelper::getBondTypeAsInt(bi);
            std::string nbr_atom_symbol = nbr_atom->getSymbol();
            replaceUncommonWithX(nbr_atom_symbol, use_full_symbols_set);

            std::string key = std::to_string(bond_type) + nbr_atom_symbol;
            if(curr_distance < dict.size())
                dict[curr_distance][key] += 1;
        }
    }
}

// Method to get atom visited order via BFS
void FingerPrintFeature::getAtomVisitOrderBFS(const RootedROMol *roMolPtr, std::vector<unsigned int> &visit_order,
                                              std::vector<unsigned int> &visit_atom_distance, int num_atoms, int depth,
                                              bool use_full_symbols_set) const {

    auto mol = roMolPtr->mol;
    auto root = roMolPtr->root;
    // get distance to root
    std::map<unsigned int, unsigned int> distances_map;
    getAtomDistanceToRoot(roMolPtr, distances_map);

    // get labels for each atoms
    std::map<unsigned int, std::string> sorting_labels;
    getSortingLabel(mol, root, distances_map, sorting_labels, use_full_symbols_set);

    // maybe a struct is a better idea
    // but this is a one off
    std::queue<const RDKit::Atom *> atom_queue;
    atom_queue.push(root);

    std::unordered_set<unsigned int> visited;
    while (!atom_queue.empty()) {
        const RDKit::Atom *curr = atom_queue.front();
        atom_queue.pop();

        // if I have see this before
        if (std::find(visit_order.begin(), visit_order.end(), curr->getIdx()) !=
            visit_order.end()) {
            continue;
        }

        auto distance_to_root = distances_map[curr->getIdx()];
        if (visit_order.size() < num_atoms && depth > distance_to_root) {
            visit_order.push_back(curr->getIdx());
            visit_atom_distance.push_back(distance_to_root);
        } else {
            break;
        }

        // use multimap since we can have duplicated labels
        std::multimap<std::string, const RDKit::Atom *> child_visit_order;
        for (auto itp = mol->getAtomNeighbors(curr); itp.first != itp.second; ++itp.first) {
            RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
            // if we have not visit this node before
            // and this node is in the visit list
            if (nbr_atom != curr) {
                std::string sorting_key = sorting_labels[nbr_atom->getIdx()];
                child_visit_order.insert(
                        std::pair<std::string, RDKit::Atom *>(sorting_key, nbr_atom));
            }
        }

        for (auto child : child_visit_order)
            atom_queue.push(child.second);
    }
}

void
FingerPrintFeature::addAdjacentMatrixRepresentation(std::vector<int> &tmp_fv, const RootedROMol *roMolPtr,
                                                    unsigned int num_atom, unsigned int depth,
                                                    bool include_adjacency_matrix, bool use_full_symbols_set) const {
    // Get visit order
    std::vector<unsigned int> visit_order;
    std::vector<unsigned int> distance;
    getAtomVisitOrderBFS(roMolPtr, visit_order, distance, num_atom, depth, false);

    if (include_adjacency_matrix)
        addAdjMatrixFeatures(tmp_fv, roMolPtr, num_atom, visit_order, distance, 0, depth, false);

    // fv.writeDebugInfo();
    // add atoms information into FP
    addAtomTypeSeqFeatures(tmp_fv, roMolPtr, num_atom, visit_order, distance, 0, depth, use_full_symbols_set);
    addDegreeFeatures(tmp_fv, roMolPtr, num_atom, visit_order);

}

void FingerPrintFeature::updateBondAtomPairDict(const RootedROMol *rootedMol, const RDKit::Atom *root,
                                                std::map<std::string, int> &dict, bool use_full_symbols_set) const {

    auto mol = rootedMol->mol;
    for (auto itp = mol->getAtomNeighbors(root); itp.first != itp.second; ++itp.first) {
        RDKit::Atom *nbr_atom = mol->getAtomWithIdx(*itp.first);
        auto bi = mol.get()->getBondBetweenAtoms(root->getIdx(), nbr_atom->getIdx());
        int bond_type = FeatureHelper::getBondTypeAsInt(bi);
        //mol.get()->getBondBetweenAtoms(root->getIdx(), nbr_atom->getIdx())->getProp("OrigBondType", bond_type);
        std::string nbr_atom_symbol = nbr_atom->getSymbol();
        replaceUncommonWithX(nbr_atom_symbol, use_full_symbols_set);

        std::string key = std::to_string(bond_type) + nbr_atom_symbol;
        dict[key] += 1;
    }
}

void FingerPrintFeature::addBondAtomPairToFeatures(std::vector<int> &tmp_fv, std::map<std::string, int> &dict,
                                                   bool no_count) const{
    for (auto const& record : dict){
        if(!no_count){
            for(int i = 0 ; i < 3; ++i){
                if(record.second > i)
                    tmp_fv.push_back(1);
                else
                    tmp_fv.push_back(0);
            }
        } else{
            if(record.second > 0)
                tmp_fv.push_back(1);
            else
                tmp_fv.push_back(0);
        }
    }
}

void FingerPrintFeature::addGenernalizedRepresentation(std::vector<int> &tmp_fv, const RootedROMol *roMolPtr,
                                                       unsigned int max_distance) const {


    auto mol = roMolPtr->mol;
    bool use_full_symbols_set = false;

    auto symbols = OKSymbolsLess();
    const int num_bond_type = 7;
    std::vector<std::map<std::string,int>> dicts(max_distance);
    //init dict
    for(auto & symbol : symbols){
        for(int bond_type = 1; bond_type <= num_bond_type; ++ bond_type){
            std::string key = std::to_string(bond_type) + symbol;
            for(int i = 0 ; i < max_distance; ++i)
                dicts[i][key] = 0;
        }
    }

    getBondAtomPairAtEachDistance(roMolPtr, dicts, use_full_symbols_set);

    int offset = 1;
    // 142 bits for atoms next to root
    for(int i = 0 ; i < offset ; ++i)
        addBondAtomPairToFeatures(tmp_fv, dicts[i], use_full_symbols_set);
    //addBondAtomPairToFeatures(tmp_fv, dicts[1], false);
    for(int i = offset; i < max_distance; ++i)
        addBondAtomPairToFeatures(tmp_fv, dicts[i], use_full_symbols_set);
}

void
FingerPrintFeature::addAdjMatrixFeatures(std::vector<int> &tmp_fv, const RootedROMol *mol, unsigned int num_atom,
                                         std::vector<unsigned int> &visit_order, std::vector<unsigned int> &distance,
                                         int min_distance, int max_distance, bool no_bond_type) const {// init a 2D vector to store matrix

    std::vector<std::vector<int>> adjacency_matrix(num_atom,
                                                   std::vector<int>(num_atom, 0));
    getAdjMatrix(mol, num_atom, visit_order, adjacency_matrix, distance, min_distance, max_distance);

    if(!no_bond_type){
        // first bit indicate if there is a bond
        // rest 7 for each bond_type, one hot encoding
        // 0, 1, ,2, 3, aromtaic as 4, Conjugated as 5
        const unsigned int num_bits_per_bond = 6;
        const unsigned int bond_type_max = num_bits_per_bond - 1;

        // we only need half of matrix exclude diagonal
        // so i start from 1 not 0
        for (int i = 0; i < num_atom; ++i) {
            for (int j = i + 1; j < num_atom; ++j) {
                std::vector<int> temp_feature(num_bits_per_bond, 0);
                // check if bond exits/defined
                if (adjacency_matrix[i][j] > 0) {
                    // one hot encoding bond type
                    int bond_type = adjacency_matrix[i][j] > bond_type_max ? bond_type_max : adjacency_matrix[i][j];
                    temp_feature[bond_type] = 1;
                }
                tmp_fv.insert(tmp_fv.end(), temp_feature.begin(), temp_feature.end());
            }
        }
    }
    else{
        for (int i = 0; i < num_atom; ++i) {
            for (int j = i + 1; j < num_atom; ++j) {
                // check if bond exits/defined
                if (adjacency_matrix[i][j] > 0)
                    tmp_fv.push_back(1);
                else
                    tmp_fv.push_back(0);
            }
        }
    }
}

void FingerPrintFeature::getAdjMatrix(const RootedROMol *mol, unsigned int num_atom,
                                      const std::vector<unsigned int> &visit_order,
                                      std::vector<std::vector<int>> &adjacency_matrix,
                                      std::vector<unsigned int> &distance,
                                      int min_distance, int max_distance) const {// make sure we only get num_atom amount of atoms, this is extra check,

    // first check is done in the getAtomVisitOrderBFS
    std::map<unsigned int, int> visit_order_map;

    for (int i = 0; i < visit_order.size() && i < num_atom; ++i){
        if((distance[i] >= min_distance)&&(distance[i] <= max_distance))
            visit_order_map[visit_order[i]] = i;
    }

    // add bound type
    for (auto bi = mol->mol->beginBonds(); bi != mol->mol->endBonds(); ++bi) {
        // for each bond find two atoms
        unsigned int begin_idx = (*bi)->getBeginAtomIdx();
        unsigned int end_idx = (*bi)->getEndAtomIdx();
        int bond_type = FeatureHelper::getBondTypeAsInt(*bi);

        // if atoms in the list
        if (visit_order_map.find(begin_idx) != visit_order_map.end() &&
            visit_order_map.find(end_idx) != visit_order_map.end()) {
            adjacency_matrix[visit_order_map[begin_idx]][visit_order_map[end_idx]] =
                    bond_type;
            adjacency_matrix[visit_order_map[end_idx]][visit_order_map[begin_idx]] =
                    bond_type;
        }
    }
}

void FingerPrintFeature::addAtomTypeSeqFeatures(std::vector<int> &tmp_fv, const RootedROMol *mol, unsigned int num_atom,
                                                const std::vector<unsigned int> &visit_order,
                                                std::vector<unsigned int> &distance, int min_distance, int max_distance,
                                                bool use_full_okay_symbol_set) const {

    unsigned int num_atom_types = use_full_okay_symbol_set ? OKsymbols().size() : OKSymbolsLess().size();

    for (int i = 0; i < num_atom; ++i) {
        std::vector<int> atom_type_feature(num_atom_types, 0);
        if (i < visit_order.size()) {
            if((distance[i] >= min_distance)&&(distance[i] <= max_distance)){
                int atom_idx = visit_order[i];
                // add atom types
                std::string symbol = mol->mol->getAtomWithIdx(atom_idx)->getSymbol();
                replaceUncommonWithX(symbol, use_full_okay_symbol_set);
                int atom_feature = getSymbolsIndex(symbol, use_full_okay_symbol_set);
                atom_type_feature[atom_feature] = 1;
            }
        }
        tmp_fv.insert(tmp_fv.end(), atom_type_feature.begin(), atom_type_feature.end());
    }
}

void FingerPrintFeature::addAtomTypeFeatures(std::vector<int> &tmp_fv, const RootedROMol *mol, unsigned int num_atom,
                                             const std::vector<unsigned int> &visit_order,
                                             const std::vector<unsigned int> &distance, bool use_full_symbols_set) const {
    const unsigned int num_atom_types = 6;
    std::vector<int> atom_type_feature(num_atom_types * 3, 0);
    for (int i = 0; i < num_atom; ++i) {
        if (i < visit_order.size()) {
            int atom_idx = visit_order[i];
            // add atom types
            std::string symbol = mol->mol->getAtomWithIdx(atom_idx)->getSymbol();
            replaceUncommonWithX(symbol, use_full_symbols_set);
            int atom_feature = getSymbolsIndex(symbol, use_full_symbols_set);
            atom_type_feature[atom_feature + distance[i] * num_atom_types] = 1;
        }
    }
    tmp_fv.insert(tmp_fv.end(), atom_type_feature.begin(), atom_type_feature.end());
}

void FingerPrintFeature::addDegreeFeatures(std::vector<int> &tmp_fv, const RootedROMol *mol, unsigned int num_atom,
                                           const std::vector<unsigned int> &visit_order) const {
    const unsigned int num_max_degree = 4;
    const unsigned int num_degree_feature_size = num_max_degree + 1;
    for (int i = 0; i < num_atom; ++i) {
        std::vector<int> atom_degree_feature(num_degree_feature_size, 0);

        if (i < visit_order.size()) {
            int atom_idx = visit_order[i];
            // add first order degree info
            int degree = mol->mol->getAtomWithIdx(atom_idx)->getDegree();
            degree = degree > num_max_degree ? num_max_degree : degree;
            atom_degree_feature[degree] = 1;
        }
        tmp_fv.insert(tmp_fv.end(), atom_degree_feature.begin(), atom_degree_feature.end());
    }
}

void FingerPrintFeature::addDistanceFeature(std::vector<int> &tmp_fv, unsigned int num_atom,
                                           const std::vector<unsigned int> &distance) const {
    const unsigned int max_distance = 10;
    const unsigned int max_features = max_distance + 1;
    for (int i = 0; i < num_atom; ++i) {
        std::vector<int> features(max_features, 0);

        if (i < distance.size()) {
            auto dis = distance[i] > max_distance ? max_distance : distance[i];
            features[dis] = 1;
        }
        tmp_fv.insert(tmp_fv.end(), features.begin(), features.end());
    }
}

// for all the samples we have max atoms with a 3 atom group is 10
// for all the samples we have max atoms with a 5 atom group is 16
// therefore  we need 50 features for arcs
void FingerPrintFeature::addAdjacentMatrixRepresentationFeature(FeatureVector &fv, const RootedROMol *mol,
                                                                unsigned int num_atom,
                                                                unsigned int max_distance,
                                                                bool include_adjacency_matrix,
                                                                bool use_full_symbols_set) const {

    std::vector<int> local_tmp_fv;
    addAdjacentMatrixRepresentation(local_tmp_fv, mol, num_atom, max_distance,
                                    include_adjacency_matrix, use_full_symbols_set);
    fv.addFeatures(local_tmp_fv);
}

void FingerPrintFeature::addGenernalizedRepresentationFeature(FeatureVector &fv, const RootedROMol *mol,
                                                              unsigned int num_atom, unsigned int max_distance) const {

    std::vector<int> local_tmp_fv;
    addGenernalizedRepresentation(local_tmp_fv, mol, max_distance);
    fv.addFeatures(local_tmp_fv);
}
