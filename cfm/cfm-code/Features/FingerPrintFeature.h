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

class FingerPrintFeature : public BreakFeature {
protected:
    // function to get part of Mol from given root
    // and all atom within given range
    // trivarse tree using BFS
    void getRemoveAtomIdxByDisatnce(romol_ptr_t mol, const RDKit::Atom *root,
                                    std::vector<unsigned int> &remove_atom_ids,
                                    int distance) const;

    void getRemoveAtomIdxByCount(romol_ptr_t mol, const RDKit::Atom *root,
                                 std::vector<unsigned int> &remove_atom_ids,
                                 int count) const;

    // remove atoms from mol in given list
    void removeAtomInTheList(RDKit::RWMol &mol,
                             std::vector<unsigned int> &remove_atom_ids) const;

    void addRDKitFingerPrintFeatures(FeatureVector &fv, const RootedROMolPtr *mol, unsigned int finger_print_size,
                                         unsigned int limitation_param, bool limited_by_distance,
                                         unsigned int finger_print_min_path, unsigned int finger_print_max_path) const;

    void addMorganFingerPrintFeatures(FeatureVector &fv, const RootedROMolPtr *mol, unsigned int finger_print_size,
                                          unsigned int path_range, int radius) const;

    void addAdjacentMatrixRepresentationFeature(FeatureVector &fv, const RootedROMolPtr *mol,
                                                unsigned int num_atom, unsigned int max_distance,
                                                bool include_adjacency_matrix) const;

    void addMorganFingerPrintFeatures(FeatureVector &fv, const RootedROMolPtr *mol,
                                      unsigned int finger_print_size, int radius) const;

    void addGenernalizedRepresentationFeature(FeatureVector &fv, const RootedROMolPtr *mol,
                                                  unsigned int num_atom, unsigned int max_distance) const;

private:

    void getAtomVisitOrderBFS(const RootedROMolPtr *roMolPtr, std::vector<unsigned int> &visit_order,
                                  std::vector<unsigned int> &visit_atom_distance, int num_atoms, int depth) const;

    std::string getSortingLabel(const romol_ptr_t mol, const RDKit::Atom *atom, const RDKit::Atom *parent_atom,
                                    std::map<unsigned int, unsigned int> &distances,
                                    std::map<unsigned int, std::string> &labels) const;

    void getAtomDistanceToRoot(const RootedROMolPtr *roMolPtr, std::map<unsigned int, unsigned int> &distances) const;

    void addMorganFingerPrint(std::vector<int> &tmp_fv, const RootedROMolPtr *mol,
                              const RDKit::Atom *root,
                              const unsigned int max_nbr_distance,
                              const unsigned int finger_print_size,
                              const int radius) const;

    void addRDKitFingerPrint(std::vector<int> &tmp_fv, const RootedROMolPtr *mol, const RDKit::Atom *root,
                             unsigned int finger_print_size, unsigned int limitation_param,
                             unsigned int finger_print_min_path, unsigned int finger_print_max_path,
                             bool limited_by_distance) const;

    void addAdjacentMatrixRepresentation(std::vector<int> &tmp_fv, const RootedROMolPtr *roMolPtr,
                                         unsigned int num_atom,
                                         unsigned int depth, bool include_con_matrix) const;

    void addGenernalizedRepresentation(std::vector<int> &tmp_fv, const RootedROMolPtr *roMolPtr,
                                           unsigned int num_atom, unsigned int max_distance) const;

    void addDegreeFeatures(std::vector<int> &tmp_fv, const RootedROMolPtr *mol, unsigned int num_atom,
                           const std::vector<unsigned int> &visit_order) const;

    void addAtomTypeSeqFeatures(std::vector<int> &tmp_fv, const RootedROMolPtr *mol, unsigned int num_atom,
                                const std::vector<unsigned int> &visit_order) const;

    void addAtomTypeFeatures(std::vector<int> &tmp_fv, const RootedROMolPtr *mol, unsigned int num_atom,
                                 const std::vector<unsigned int> &visit_order, const std::vector<unsigned int> &distance,
                                unsigned int max_distance) const;

    void addAdjMatrixFeatures(std::vector<int> &tmp_fv, const RootedROMolPtr *mol, unsigned int num_atom,
                              std::vector<unsigned int> &visit_order) const;

    void addDistanceFeature(std::vector<int> &tmp_fv, unsigned int num_atom, const std::vector<unsigned int> &distance) const;
};
