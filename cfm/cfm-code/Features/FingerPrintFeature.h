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
    void getRemoveAtomIdxByDisatnce(romol_ptr_t mol, const RDKit::Atom *root,
                                    std::vector<unsigned int> &remove_atom_ids,
                                    int distance) const;

    void getRemoveAtomIdxByCount(romol_ptr_t mol, const RDKit::Atom *root,
                                    std::vector<unsigned int> &remove_atom_ids,
                                    int count) const;

    // remove atoms from mol in given list
    void removeAtomInTheList(RDKit::RWMol &mol,
                             std::vector<unsigned int> &remove_atom_ids) const;

    // replace bond type with orig bond type
    void replaceWithOrigBondType(RDKit::RWMol &mol) const;

    void addRDKitFingerPrintFeatures(
            FeatureVector &fv, const RootedROMolPtr *mol,
            unsigned int finger_print_size, unsigned int max_nbr_distance,
            int ring_break, unsigned int finger_print_min_path,
            unsigned int finger_print_max_path) const;

    void addMorganFingerPrintFeatures(FeatureVector &fv,
                                      const RootedROMolPtr *mol,
                                      unsigned int finger_print_size,
                                      unsigned int path_range,
                                      int ring_break,
                                      int radius) const;

    void addAdjacentMatrixRepresentationFeature(FeatureVector &fv, const RootedROMolPtr *mol,
                                                unsigned int num_atom, int ring_break,
                                                bool include_adjacency_matrix) const;

private:
    // void getAtomsWithRange(int range);

    // void getAtomVisitOrderDFS(const romol_ptr_t mol, const RDKit::Atom *atom,
    //                          const RDKit::Atom *prev_atom, int range,
    //                          std::vector<unsigned int> &visited) const;
    enum Bfs_Stop_Logic {
        RANGE_ONLY = 0, SIZE_ONLY, RANGE_OR_SIZE
    };

    void getAtomVisitOrderBFS(const romol_ptr_t mol, const RDKit::Atom *root,
                              std::vector<unsigned int> &visit_order, int num_atoms,
                              std::map<int, int> &path_record) const;

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
                              unsigned int max_nbr_distance,
                              unsigned int finger_print_size,
                              int radius) const;

    void addRDKitFingerPrint(FeatureVector &fv, const RootedROMolPtr *mol,
                             const RDKit::Atom *root,
                             unsigned int finger_print_size,
                             unsigned int max_nbr_distance,
                             unsigned int finger_print_min_path,
                             unsigned int finger_print_max_path) const;

    void addAdjacentMatrixRepresentation(FeatureVector &fv, const RootedROMolPtr *mol,
                                         const RDKit::Atom *root, unsigned int num_atom,
                                         bool include_con_matrix) const;
};
