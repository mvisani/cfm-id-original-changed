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

#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/MolOps.h>

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
 
  
  void FingerPrintFeature::getAtomVisitOrder(
    const romol_ptr_t mol, const RDKit::Atom *atom,
    const RDKit::Atom *prev_atom, std::vector<unsigned int> &visited) const {

    if (atom == nullptr ) {
      return;
    }
    if (visited.find(atom->getIdx()) != visited.end()) {
      return;
    }

    if (range == 0) {
       return;
    }
        
    visited.push_back(atom->getIdx());

    // get all the neighbors
    std::vector<RDKit::Atom *> childs;
    for (auto itp = mol.get()->getAtomNeighbors(atom); itp.first != itp.second;
        ++itp.first) {
      RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
      if (nbr_atom != prev_atom) {
        childs.push_back(nbr_atom);
      }
    }

    // sort child, we need a visiting order
    std::sort(childs.begin(), childs.end(), 
    [](RDKit::Atom * a, RDKit::Atom * b) -> bool
    { 
        if(a->getAtomicNum() ==  b->getAtomicNum())
        {
          return a->getDegree() > b->getDegree();
        }
        else
        {
          return a->getAtomicNum() > b->getAtomicNum();
        }
      });
    
    for(auto itp = childs.begin(); itp != childs.end(); ++itp )
    {
      getAtomVisitOrder(mol, itp, atom, visited); 
    }
  }
 
  /*
  flatten: ajcent matrix[]
  first bit - if there is a bond
  follwing 4 bits hot for type

  list of atom per node

  list of degrees per node

  */

  void FingerPrintFeature::addHomeMadeFingerPrint(std::vector<int> &ids, 
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
    std::vector<unsigned int> visitOrder;
    getAtomVisitOrder(mol->mol, mol->root, nullptr, visitOrder);
    
    // for all the bounds
    for (auto bi = mol->mol.get().beginBonds(); bi != mol->mol.get().endBonds(); ++bi) {
      // for each bond
      
    }
    
    // add features
    for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
      fv.addFeature((*fingerPrint)[i]);
    }
  
    if (ring_break) {
    
    } else {
      for (unsigned int i = 0; i < finger_print_size; ++i) {
        fv.addFeature(0);
      }
    }
}