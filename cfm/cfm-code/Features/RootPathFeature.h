/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# RootPathFeature.h
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
#include <GraphMol/RWMol.h>

class RootPathFeature : public Feature {
protected:
  typedef std::vector<std::string> path_t;
  // function to compute path with given length from a root
  void computeRootPaths(std::vector<path_t> &paths, const RootedROMolPtr *mol,
                        int len, bool ring_break, bool with_bond) const;

  // function to add features with a length of two
  void addRootPairFeatures(FeatureVector &fv, std::vector<path_t> &paths,
                           int ring_break) const;

  // function to add features with a length of three
  void addRootTripleFeatures(FeatureVector &fv, std::vector<path_t> &paths,
                             int ring_break, bool with_bond) const;

  // function to get part of Mol from given root
  // and all atom within given range
  void getRmoveAtomIdxOfRange(const romol_ptr_t mol, const RDKit::Atom *atom,
                              const RDKit::Atom *prev_atom,
                              std::vector<unsigned int> &remove_atom_ids,
                              int range) const;

  void removeAtomNotInTheList(
      RDKit::RWMol &mol,
      const std::vector<unsigned int> &remove_atom_ids) const;

private:
  // function to add path from given atom
  void addPathsFromAtom(std::vector<path_t> &paths, const RDKit::Atom *atom,
                        const romol_ptr_t mol, const RDKit::Atom *prev_atom,
                        path_t &path_so_far, int len, bool with_bond) const;
};
