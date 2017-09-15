/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# RingFeatures.cpp
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
#include "RingFeatures.h"

void RingFeatures::compute(FeatureVector &fv, const RootedROMolPtr *ion,
                           const RootedROMolPtr *nl) const {

  // Non-Ring break
  int ring_break;
  nl->mol.get()->getProp("IsRingBreak", ring_break);
  if (!ring_break) {
    for (int i = 0; i < 12; i++)
      fv.addFeature(0);
  } else {
    // Ring Break
    // Find a broken bond, and fetch the aromaticity labels from it
    int is_arom, is_dbl_arom;
    nl->mol.get()->getProp("IsAromaticRingBreak", is_arom);
    nl->mol.get()->getProp("IsAromaticDblRingBreak", is_dbl_arom);

    int root_dist_ion = calcRootDistance(ion);
    int root_dist_nl = calcRootDistance(nl);

    int break_dist = std::min(root_dist_ion, root_dist_nl) + 1;
    int ring_size = root_dist_nl + root_dist_ion + 2;

    fv.addFeature(!is_arom);
    fv.addFeature(is_arom);
    fv.addFeature(is_dbl_arom);
    for (int dist = 1; dist <= 3; dist++)
      fv.addFeature(break_dist == dist);
    fv.addFeature(break_dist > 3);
    for (int rsize = 3; rsize <= 6; rsize++)
      fv.addFeature(rsize == ring_size);
    fv.addFeature(ring_size > 6);
  }
}

int RingFeatures::calcRootDistance(const RootedROMolPtr *mol) const {
  // Compute the distance from the root to the other root (breadth-first-search)
  typedef std::pair<RDKit::Atom *, int> item_t;
  std::set<int> seen_idxs;
  std::queue<item_t> queue;
  queue.push(item_t(mol->root, 0));
  while (queue.size() > 0) {
    item_t item = queue.front();
    RDKit::Atom *at = item.first;
    int depth = item.second;
    seen_idxs.insert(at->getIdx());

    // Have we found the root?
    if (at == mol->other_root)
      return depth;

    // If not, keep looking
    RDKit::ROMol::ADJ_ITER_PAIR itp = mol->mol.get()->getAtomNeighbors(at);
    for (; itp.first != itp.second; ++itp.first) {
      RDKit::Atom *nbr_atom = mol->mol.get()->getAtomWithIdx(*itp.first);
      if (seen_idxs.find(nbr_atom->getIdx()) == seen_idxs.end())
        queue.push(item_t(nbr_atom, depth + 1));
    }
    queue.pop();
  }
  return -1; // OtherRoot not found
}