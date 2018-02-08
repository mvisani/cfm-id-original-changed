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
#pragma once

#include "../Feature.h"

#include <unordered_set>

#include <GraphMol/RWMol.h>

class FingerPrintFeature : public Feature {
protected:
  // function to get part of Mol from given root
  // and all atom within given range
  // trivarse tree using BFS
  void getRemoveAtomIdx(const romol_ptr_t mol, const RDKit::Atom *root,
                        std::vector<unsigned int> &remove_atom_ids,
                        int distance) const;

  // remove atoms from mol in given list
  void removeAtomInTheList(RDKit::RWMol &mol,
                           std::vector<unsigned int> &remove_atom_ids) const;

  // replace bond type with orig bond type
  void replaceWithOrigBondType(RDKit::RWMol &mol) const;

  void addRDKitFingerPrintFeatures(
      FeatureVector &fv, const RootedROMolPtr *mol,
      const unsigned int finger_print_size, const unsigned int max_nbr_distance,
      const int ring_break, const unsigned int finger_print_min_path,
      const unsigned int finger_print_max_path) const;

  void addMorganFingerPrintFeatures(FeatureVector &fv,
                                    const RootedROMolPtr *mol,
                                    const unsigned int finger_print_size,
                                    const unsigned int path_range,
                                    const int ring_break,
                                    const int radius) const;

  void addAdjacentMatrixRepresentationFeature(
          FeatureVector &fv, const RootedROMolPtr *mol,
          const unsigned int max_nbr_distance, const unsigned int num_atom,
          const int ring_break) const;

private:
  // void getAtomsWithRange(int range);

  // void getAtomVisitOrderDFS(const romol_ptr_t mol, const RDKit::Atom *atom,
  //                          const RDKit::Atom *prev_atom, int range,
  //                          std::vector<unsigned int> &visited) const;
  enum Bfs_Stop_Logic { RANGE_ONLY = 0, SIZE_ONLY, RANGE_OR_SIZE };

  void getAtomVisitOrderBFS(const romol_ptr_t mol, const RDKit::Atom *root,
                            std::vector<unsigned int> &visit_order, int range,
                            int size, Bfs_Stop_Logic stop_logic) const;

  std::string getSortingLabel(const romol_ptr_t mol, const RDKit::Atom *atom,
                              const RDKit::Atom *parent_atom, int depth) const;

  /*
  std::string
  getSortingLabels(const romol_ptr_t mol, const RDKit::Atom *atom,
                   const RDKit::Atom *prev_atom, int range,
                   std::unordered_set<unsigned int> &visited,
                   std::map<unsigned int, std::string> &sorting_labels) const;
  */
  /*void addMorganFingerPrint(FeatureVector &fv,
                            const RootedROMolPtr *mol,
                            const RDKit::Atom *root,
                            const unsigned int path_range,
                            const int radius) const;*/

  void addMorganFingerPrint(FeatureVector &fv, const RootedROMolPtr *mol,
                            const RDKit::Atom *root,
                            const unsigned int max_nbr_distance,
                            const unsigned int finger_print_size,
                            const int radius) const;

  void addRDKitFingerPrint(FeatureVector &fv, const RootedROMolPtr *mol,
                           const RDKit::Atom *root,
                           const unsigned int finger_print_size,
                           const unsigned int max_nbr_distance,
                           const unsigned int finger_print_min_path,
                           const unsigned int finger_print_max_path) const;

  void addAdjacentMatrixRepresentation(FeatureVector &fv,
                                       const RootedROMolPtr *mol,
                                       const RDKit::Atom *root,
                                       const unsigned int max_nbr_distance,
                                       const unsigned int num_atom) const;
};
