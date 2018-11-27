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

    void addAdjacentMatrixRepresentationFeature(FeatureVector &fv, const RootedROMolPtr *mol, unsigned int num_atom,
                                                    bool include_adjacency_matrix) const;

    void addMorganFingerPrintFeatures(FeatureVector &fv, const RootedROMolPtr *mol,
                                      unsigned int finger_print_size, int radius) const;

private:

    void getAtomVisitOrderBFS(const romol_ptr_t mol, const RDKit::Atom *root,
                                  std::vector<unsigned int> &visit_order, int num_atoms) const;

    std::string getSortingLabel(const romol_ptr_t mol, const RDKit::Atom *atom,
                                const RDKit::Atom *parent_atom, int depth) const;


    void addMorganFingerPrint(std::vector<int> &tmp_fv, const RootedROMolPtr *mol,
                              const RDKit::Atom *root,
                              const unsigned int max_nbr_distance,
                              const unsigned int finger_print_size,
                              const int radius) const;

    void addRDKitFingerPrint(std::vector<int> &tmp_fv, const RootedROMolPtr *mol, const RDKit::Atom *root,
                             unsigned int finger_print_size, unsigned int limitation_param,
                             unsigned int finger_print_min_path, unsigned int finger_print_max_path,
                             bool limited_by_distance) const;

    void addAdjacentMatrixRepresentation(std::vector<int> &tmp_fv, const RootedROMolPtr *mol,
                                         const RDKit::Atom *root, unsigned int num_atom,
                                         bool include_con_matrix) const;

};
