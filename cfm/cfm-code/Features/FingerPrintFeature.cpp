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

#include <map>

void FingerPrintFeature::getRemoveAtomIdxOfRange(
    const romol_ptr_t mol, const RDKit::Atom *atom,
    const RDKit::Atom *prev_atom, std::vector<unsigned int> &remove_atom_ids,
    std::unordered_set<unsigned int> &visited, int range) const {

  if (atom == nullptr) {
    return;
  }
  
  if (visited.find(atom->getIdx()) != visited.end()) {
    return;
  }

  visited.insert(atom->getIdx());

  // if range <= 0 , that means this atom is out of range
  // thus we need remove it
  // also we does not care Hs
  if (range <= 0 ||
      atom->getSymbol() == "H") 
  {
    remove_atom_ids.push_back(atom->getIdx());
  }

  // get all the neighbors
  for (auto itp = mol.get()->getAtomNeighbors(atom); itp.first != itp.second;
       ++itp.first) {
    RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
    if (nbr_atom != prev_atom) {
      getRemoveAtomIdxOfRange(mol, nbr_atom, atom, remove_atom_ids, visited,
                              range - 1);
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

void FingerPrintFeature::addRDKitFingerPrint(FeatureVector &fv, const RootedROMolPtr *mol,
                                            const unsigned int finger_print_size,
                                            const unsigned int path_range, const int ring_break,
                                            const unsigned int finger_print_min_path,
                                            const unsigned int finger_print_max_path) const {

    RDKit::ROMol &nl_ref = *(mol->mol.get());
    // Get list of atom we need to remove
    std::vector<unsigned int> remove_atom_ids;
    std::unordered_set<unsigned int> visited;
    getRemoveAtomIdxOfRange(mol->mol, mol->root, nullptr, remove_atom_ids,
                            visited, path_range);
  
    // Get Mol Object and remove atoms
    RDKit::RWMol part;
    part.insertMol(*(mol->mol.get()));
    removeAtomInTheList(part, remove_atom_ids);
    RDKit::MolOps::sanitizeMol(part);
  
    // replace bond with OrigBondType
    // replaceWithOrigBondType(part);
  
    // Get finger prints with size
    ExplicitBitVect *fingerPrint = RDKit::RDKFingerprintMol(
        part, finger_print_min_path, finger_print_max_path, finger_print_size);
  
    for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
      fv.addFeature((*fingerPrint)[i]);
    }
  
    if (ring_break) {
      remove_atom_ids.clear();
      getRemoveAtomIdxOfRange(mol->mol, mol->other_root, nullptr, remove_atom_ids,
                              visited, path_range);
  
      // Get Mol Object and remove atoms
      RDKit::RWMol other_part;
      other_part.insertMol(*(mol->mol.get()));
      this->removeAtomInTheList(other_part, remove_atom_ids);
  
      // we don't want to part to be santilized
      // because if we are taking part of ring
      RDKit::MolOps::sanitizeMol(part);
  
      ExplicitBitVect *fingerPrint =
      RDKit::RDKFingerprintMol(other_part, finger_print_min_path,
                                   finger_print_max_path, finger_print_size);
  
      for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
        fv.addFeature((*fingerPrint)[i]);
      }
    } else {
      for (unsigned int i = 0; i < finger_print_size; ++i) {
        fv.addFeature(0);
      }
    }
}

void FingerPrintFeature::addMorganFingerPrint(FeatureVector &fv, 
                                              const RootedROMolPtr *mol,
                                              const unsigned int finger_print_size,
                                              const unsigned int path_range, 
                                              const int ring_break,
                                              const int radius) const {
    
    RDKit::ROMol &nl_ref = *(mol->mol.get());
    // Get list of atom we need to remove
    std::vector<unsigned int> remove_atom_ids;
    std::unordered_set<unsigned int> visited;
    getRemoveAtomIdxOfRange(mol->mol, mol->root, nullptr, remove_atom_ids,
                            visited, path_range);
  
    // Get Mol Object and remove atoms
    RDKit::RWMol part;
    part.insertMol(*(mol->mol.get()));
    removeAtomInTheList(part, remove_atom_ids);
    
    // replace bond with OrigBondType
    // replaceWithOrigBondType(part);
  
    // Get finger prints with size
    ExplicitBitVect *fingerPrint = RDKit::MorganFingerprints::getFingerprintAsBitVect (
        part, radius, finger_print_size);
  
    for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
      fv.addFeature((*fingerPrint)[i]);
    }
  
    if (ring_break) {
      remove_atom_ids.clear();
      getRemoveAtomIdxOfRange(mol->mol, mol->other_root, nullptr, remove_atom_ids,
                              visited, path_range);
  
      // Get Mol Object and remove atoms
      RDKit::RWMol other_part;
      other_part.insertMol(*(mol->mol.get()));
      this->removeAtomInTheList(other_part, remove_atom_ids);
  
      // we don't want to part to be santilized
      // because if we are taking part of ring
      // RDKit::MolOps::sanitizeMol(part);
  
      ExplicitBitVect *fingerPrint =
          RDKit::MorganFingerprints::getFingerprintAsBitVect(other_part, radius, finger_print_size);
  
      for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
        fv.addFeature((*fingerPrint)[i]);
      }
    } else {
      for (unsigned int i = 0; i < finger_print_size; ++i) {
        fv.addFeature(0);
      }
    }
  }
 
  std::string FingerPrintFeature::getSortingLabels(
    const romol_ptr_t mol, const RDKit::Atom *atom,
    const RDKit::Atom *prev_atom,  int range,
    std::unordered_set<unsigned int> & visited,
    std::map<unsigned int, std::string> &sorting_labels) const
  {
      if (atom == nullptr  
          || range == 0 
          || visited.find(atom->getIdx()) != visited.end()) {
        return "";
      }

    visited.insert(atom->getIdx());

    // get all the neighbors
    std::vector<std::string> children_keys;
    for (auto itp = mol.get()->getAtomNeighbors(atom); itp.first != itp.second;
        ++itp.first) {
      RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
      if (nbr_atom != prev_atom) {
        std::string child_key = "";
        child_key = getSortingLabels(mol, nbr_atom, atom, range-1, visited, sorting_labels); 
        children_keys.push_back(child_key);
      }
    }

    // we don't care root node
    // there is no need to figure out which root to visit first anyway
    std::string atom_key = "";
    if(prev_atom != nullptr)
    {
      std::sort(children_keys.begin(), children_keys.end());
      // add bond type
      int bond_int = FeatureHelper::getBondTypeAsInt(mol.get()->getBondBetweenAtoms(prev_atom->getIdx(),atom->getIdx()));
      atom_key += std::to_string(bond_int);

      // add atom symbol:
      // TODO: figure out which way is better 
      // replace symbol here or use true symbol
      std::string symbol_str = atom->getSymbol();
      replaceUncommonWithX(symbol_str);
      atom_key += symbol_str;

      // get child atom keys str
      std::string children_atom_key = "";
      for(auto child_key : children_keys)
      {
        children_atom_key += child_key;
      }
      sorting_labels.insert(std::pair<unsigned int, std::string>(atom->getIdx(),atom_key + children_atom_key));
    }
    return atom_key;
  }

  //Method to get atom visited order
  void FingerPrintFeature::getAtomVisitOrder(
    const romol_ptr_t mol, const RDKit::Atom *atom,
    const RDKit::Atom *prev_atom,  int range,
    std::vector<unsigned int> &visited,
    const std::map<unsigned int, std::string> & sorting_labels) const {
    
    if (atom == nullptr  
          || range == 0 
          || std::find(visited.begin(), visited.end(), atom->getIdx()) != visited.end()) {
        return;
    }
    visited.push_back(atom->getIdx());

    // create a map to visit child
    // use multimap since we can have duplicated labels
    std::multimap<std::string, RDKit::Atom *> child_visit_order;

    // get all the neighbors and insert into a multimap
    // use map to do sorting 
    std::vector<RDKit::Atom *> childs;
    for (auto itp = mol.get()->getAtomNeighbors(atom); itp.first != itp.second;
        ++itp.first) {
      RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
      if (nbr_atom != prev_atom) {
        std::string sorting_key = sorting_labels.find(nbr_atom->getIdx())->second;
        child_visit_order.insert(std::pair<std::string, RDKit::Atom *>(sorting_key, nbr_atom));
      }
    }

    for(auto child = child_visit_order.begin(); child != child_visit_order.end(); ++child )
    {
      getAtomVisitOrder(mol, child->second, atom, range-1, visited, sorting_labels); 
    }
  }
 
  /*
  flatten: ajcent matrix[]
  first bit - if there is a bond
  follwing 4 bits hot for type

  list of atom per node

  list of degrees per node

  */

  void FingerPrintFeature::addAdjacentMatrixRepesentation(FeatureVector &fv, 
                                              const RootedROMolPtr *mol,
                                              const RDKit::Atom *root,
                                              const unsigned int path_range,
                                              const unsigned int num_atom) const {
    
    RDKit::ROMol &nl_ref = *(mol->mol.get());
    
    // Get list of atom we need to remove
    //std::vector<unsigned int> remove_atom_ids;

    //getRemoveAtomIdxOfRange(mol->mol, root, nullptr, remove_atom_ids,
    //                       visited, path_range);
  
    // Get Mol Object and remove atoms
    // RDKit::RWMol part;
    // part.insertMol(*(mol->mol.get()));
    // removeAtomInTheList(part, remove_atom_ids);
    
    // Get labels for find visit order
    std::unordered_set<unsigned int> visited;
    std::map<unsigned int, std::string> sorting_labels;
    getSortingLabels(mol->mol, root, nullptr, path_range, visited, sorting_labels);

    // Get visit order
    std::vector<unsigned int> visit_order;
    getAtomVisitOrder(mol->mol, root, nullptr, path_range, visit_order, sorting_labels);
    

    std::map<unsigned int, int> visit_order_map;
    
    for(int i = 0; i < visit_order.size(); ++i )
    {
      visit_order_map[visit_order[i]] = i;
    }

    int ajcent_matrix[num_atom][num_atom] = { {0} };

    for (auto bi = mol->mol->beginBonds(); bi != mol->mol->endBonds(); ++bi) {
    // for each bond
      unsigned int beginIdx = (*bi)->getBeginAtomIdx();
      unsigned int endIdx = (*bi)->getEndAtomIdx();
      int bond_type = FeatureHelper::getBondTypeAsInt(*bi);

      // if atoms in the list
      if(visit_order_map.find(beginIdx) != visit_order_map.end() 
        && visit_order_map.find(endIdx) != visit_order_map.end())
        {
            ajcent_matrix[visit_order_map[beginIdx]][visit_order_map[endIdx]] = bond_type;
            ajcent_matrix[visit_order_map[endIdx]][visit_order_map[beginIdx]] = bond_type;
        }
    }
    
    // we only need half of matrix miuns diagonal
    // so i start from 1 not 0
    for (int i = 1 ; i < num_atom; ++i)
    {
      for(int j = i + 1; j < num_atom; ++j)
      {
          int temp_feature[5] = {0};
          if(ajcent_matrix[i][j] > 0)
          {
            // first bit indicate if there is a bond
            temp_feature[0] = 1;
            // one hot encoding bond type
            temp_feature[ajcent_matrix[i][j]] = 1;
          }

           // TODO Change to C++11 array
          fv.addFeatures(temp_feature, 5);
      }
    }

    for(int i = 0 ; i < num_atom; ++i)
    {
      int atom_type_feature[5] = {0};
      int atom_degree_feature[5] = {0};
      if(i < visit_order.size())
      {
        int atom_idx = visit_order[i];
        std::string symbol = mol->mol->getAtomWithIdx(atom_idx)->getSymbol();
        int degree = mol->mol->getAtomWithIdx(atom_idx)->getDegree();
        degree = degree > 4 ? 4 : degree;
        int atom_feature = getSymbolsLessIndex(symbol);
        atom_type_feature[atom_feature] = 1;
        atom_degree_feature[degree] = 1;
      }
      // TODO Change to C++11 array
      fv.addFeatures(atom_type_feature, 5);
      fv.addFeatures(atom_degree_feature, 5);
    }
  }


  // for all the samples we have max atoms with a 3 atom group is 10
  // for all the samples we have max atoms with a 3 atom group is 16
  // therefore  we need 50 features for arcs
  void FingerPrintFeature::addAdjacentMatrixRepesentationFeature(FeatureVector &fv, 
                                                          const RootedROMolPtr *mol,
                                                          const unsigned int path_range, 
                                                          const unsigned int num_atom,
                                                          const int ring_break) const {

    addAdjacentMatrixRepesentation(fv,mol, mol->root, path_range, num_atom);
    if(ring_break > 0)
    {
      addAdjacentMatrixRepesentation(fv, mol, mol->other_root, path_range, num_atom);
    }
    else
    {
      unsigned int feature_size = num_atom * (num_atom-1) / 2 * 5 + num_atom * 10;
      double empty_feature[feature_size] = {0};
      fv.addFeatures(empty_feature, feature_size);
    }
  }