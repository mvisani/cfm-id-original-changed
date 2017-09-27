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

#include <algorithm>


#include <DataStructs/ExplicitBitVect.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/BondIterators.h>

void RootPathFeature::computeRootPaths(std::vector<path_t> &paths,
                                       const RootedROMolPtr *mol, int len,
                                       bool ring_break,
                                       bool with_bond = false) const {
  path_t path_so_far;
  addPathsFromAtom(paths, mol->root, mol->mol, mol->root, path_so_far, len,
                   with_bond);
  if (ring_break) {
    addPathsFromAtom(paths, mol->other_root, mol->mol, mol->other_root,
                     path_so_far, len, with_bond);
  }
}

void RootPathFeature::addPathsFromAtom(std::vector<path_t> &paths,
                                       const RDKit::Atom *atom,
                                       const romol_ptr_t mol,
                                       const RDKit::Atom *prev_atom,
                                       path_t &path_so_far, int len,
                                       bool with_bond) const {
  // Add the current symbol
  std::string symbol = atom->getSymbol();
  replaceUncommonWithX(symbol);
  path_so_far.push_back(symbol);

  // Iterate until len is reached, then add the path
  if (len > 1) {
    RDKit::ROMol::ADJ_ITER_PAIR itp = mol.get()->getAtomNeighbors(atom);
    for (; itp.first != itp.second; ++itp.first) {
      RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
      if (with_bond == true) {
        int bond_type = 0;
        mol.get()
            ->getBondBetweenAtoms(atom->getIdx(), nbr_atom->getIdx())
            ->getProp("OrigBondType", bond_type);
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

void RootPathFeature::getRemoveAtomIdxOfRange(const romol_ptr_t mol, 
                                              const RDKit::Atom *atom,
                                              const RDKit::Atom *prev_atom, 
                                              std::vector<unsigned int> &remove_atom_ids,
                                              std::unordered_set<unsigned int> &visited,
                                              int range) const {

  if (atom == nullptr)
  {
    return;
  }
  if (visited.find(atom->getIdx()) != visited.end())
  {
    return;
  }

  visited.insert(atom->getIdx());
  
  // if range <= 0 , that means this atom is out of range
  // thus we need remove it
  if (range <= 0) {
    remove_atom_ids.push_back(atom->getIdx());
  }

  // get all the neighbors
  for (auto itp = mol.get()->getAtomNeighbors(atom); itp.first != itp.second;
       ++itp.first) {
    RDKit::Atom *nbr_atom = mol.get()->getAtomWithIdx(*itp.first);
    if (nbr_atom != prev_atom) {
      getRemoveAtomIdxOfRange(mol, nbr_atom, atom, remove_atom_ids,visited, range - 1);
    }
  }
}

void RootPathFeature::removeAtomNotInTheList(
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
void RootPathFeature::replaceWithOrigBondType(RDKit::RWMol &rwmol) const {

  for (auto bi = rwmol.beginBonds(); bi != rwmol.endBonds(); ++bi) {
    int bond_type;
    (*bi)->getProp("OrigBondTypeRaw", bond_type);
    (*bi)->setBondType(static_cast<RDKit::Bond::BondType>(bond_type));
  }
}

void RootPathFeature::addFingerPrint(FeatureVector &fv,
                                     const RootedROMolPtr *mol,
                                     const int finger_print_size,
                                     const int path_range,
                                     int ring_break) const {

  RDKit::ROMol &nl_ref = *(mol->mol.get());
  // Get list of atom we need to remove
  std::vector<unsigned int> remove_atom_ids;
  std::unordered_set<unsigned int> visited;
  getRemoveAtomIdxOfRange(mol->mol,
                          mol->root, 
                          nullptr, 
                          remove_atom_ids,
                          visited,
                          path_range);

  // Get Mol Object and remove atoms
  RDKit::RWMol part;
  part.insertMol(*(mol->mol.get()));
  removeAtomNotInTheList(part, remove_atom_ids);
  // replace bond with OrigBondType
  replaceWithOrigBondType(part);

  // Get finger prints with size
  unsigned int minPath = 1;
  unsigned int maxPath = 7;
  ExplicitBitVect *fingerPrint =
      RDKit::RDKFingerprintMol(part, minPath, maxPath, finger_print_size);

  for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
    fv.addFeature((*fingerPrint)[i]);
  }

  if (ring_break) {
    remove_atom_ids.clear();
    getRemoveAtomIdxOfRange(mol->mol,
                            mol->root, 
                            nullptr, 
                            remove_atom_ids,
                            visited,
                            path_range);

    // Get Mol Object and remove atoms
    RDKit::RWMol other_part;
    part.insertMol(*(mol->mol.get()));
    this->removeAtomNotInTheList(other_part, remove_atom_ids);

    // we don't want to part to be santilized
    // because if we are taking part of ring
    // RDKit::MolOps::sanitizeMol(part);

    ExplicitBitVect *fingerPrint = RDKit::RDKFingerprintMol(
        other_part, minPath, maxPath, finger_print_size);
    for (unsigned int i = 0; i < fingerPrint->getNumBits(); ++i) {
      fv.addFeature((*fingerPrint)[i]);
    }
  } else {
    for (unsigned int i = 0; i < finger_print_size; ++i) {
      fv.addFeature(0);
    }
  }
}

void RootPathFeature::addRootPairFeatures(FeatureVector &fv,
                                          std::vector<path_t> &paths,
                                          int ring_break) const {
  // Add a feature indicating that there are no pairs
  fv.addFeature((double)(paths.size() == 0));

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
  fv.addFeature((double)(paths.size() == 0));

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
            num_symbol_type * feature_idx_offset + getSymbolsLessIndex(symbol);
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
