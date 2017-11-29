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
    void getRemoveAtomIdxOfRange(const romol_ptr_t mol, const RDKit::Atom *atom,
                                 const RDKit::Atom *prev_atom,
                                 std::vector<unsigned int> &remove_atom_ids,
                                 std::unordered_set<unsigned int> &visited,
                                 int range) const;

    // remove atoms from mol in given list
    void removeAtomInTheList(RDKit::RWMol &mol,
                             std::vector<unsigned int> &remove_atom_ids) const;

    // replace bond type with orig bond type
    void replaceWithOrigBondType(RDKit::RWMol &mol) const;

    void addRDKitFingerPrint(FeatureVector &fv, const RootedROMolPtr *mol,
                             const unsigned int finger_print_size,
                             const unsigned int path_range, const int ring_break,
                             const unsigned int finger_print_min_path,
                             const unsigned int finger_print_max_path) const;

    void addMorganFingerPrint(FeatureVector &fv,
                              const RootedROMolPtr *mol,
                              const unsigned int finger_print_size,
                              const unsigned int path_range,
                              const int ring_break,
                              const int radius) const;


    void addAdjacentMatrixRepesentationFeature(FeatureVector &fv,
                                               const RootedROMolPtr *mol,
                                               const unsigned int path_range,
                                               const unsigned int num_atom,
                                               const int ring_break) const;

private:
    void getAtomsWithRange(int range);

    void getAtomVisitOrderDFS(const romol_ptr_t mol, const RDKit::Atom *atom,
                           const RDKit::Atom *prev_atom, int range, 
                           std::vector<unsigned int> &visited) const;
    
    void getAtomVisitOrderBFS(
        const romol_ptr_t mol, const RDKit::Atom *root, int range,
        std::vector<unsigned int> &visit_order) const;
    
    std::string getSortinglabel(
    const romol_ptr_t mol, const RDKit::Atom *atom,
    const RDKit::Atom *parent_atom, bool include_child = true) const;

    std::string getSortingLabels(
            const romol_ptr_t mol, const RDKit::Atom *atom,
            const RDKit::Atom *prev_atom, int range,
            std::unordered_set<unsigned int> &visited,
            std::map<unsigned int, std::string> &sorting_labels) const;

    void addAdjacentMatrixRepresentation(FeatureVector &fv,
                                         const RootedROMolPtr *mol,
                                         const RDKit::Atom *root,
                                         const unsigned int path_range,
                                         const unsigned int num_atom) const;

};
