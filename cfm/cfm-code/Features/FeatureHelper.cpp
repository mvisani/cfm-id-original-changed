/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FeatureHelper.cpp
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
#include "FeatureHelper.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/Substruct/SubstructMatch.h>

void FeatureHelper::initialiseRoots(RDKit::RWMol *rwmol) {
    RDKit::ROMol::AtomIterator ai;
    for (ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai)
        (*ai)->setProp("Root", 0);
}

void FeatureHelper::labelGasteigers(RDKit::RWMol *rwmol) {
    // For each atom, will store the result in prop "OrigGasteigerCharge"
    RDKit::computeGasteigerCharges(rwmol);
    RDKit::ROMol::AtomIterator ai;
    for (ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai) {
        double gc;
        (*ai)->getProp<double>("_GasteigerCharge", gc);
        (*ai)->setProp("OrigGasteigerCharge", gc);
    }
}

void FeatureHelper::labelFunctionalGroups(RDKit::RWMol *rwmol, bool extra) {

    const RDKit::MOL_SPTR_VECT &fgrps = fparams->getFuncGroups();
    const RDKit::MOL_SPTR_VECT &xfgrps = xfparams->getFuncGroups();

    std::vector<std::vector<unsigned int>> atom_fg_idxs(rwmol->getNumAtoms());
    std::vector<std::string> atom_fg_strs(rwmol->getNumAtoms(), "");

    std::string prop_name;
    int num_grps, idx = 0;
    RDKit::MOL_SPTR_VECT::const_iterator fgrpi, fgrpe;
    if (extra) {
        fgrpi = xfgrps.begin();
        fgrpe = xfgrps.end();
        num_grps = NUM_EXTRA_FGRPS;
        prop_name = "ExtraFunctionalGroups";
    } else {
        fgrpi = fgrps.begin();
        fgrpe = fgrps.end();
        num_grps = NUM_FGRPS;
        prop_name = "FunctionalGroups";
    }

    for (; fgrpi != fgrpe; ++fgrpi, idx++) {
        std::string fname;
        (*fgrpi)->getProp("_Name", fname);
        std::vector<RDKit::MatchVectType>
                fgpMatches; // The format for each match is (queryAtomIdx, molAtomIdx)
        RDKit::SubstructMatch(*rwmol, *(fgrpi->get()), fgpMatches);

        std::vector<RDKit::MatchVectType>::const_iterator mat_it =
                fgpMatches.begin();
        for (; mat_it != fgpMatches.end(); ++mat_it){
            for (auto it = (*mat_it).begin(); it != (*mat_it).end(); ++it){
                atom_fg_idxs[it->second].push_back(idx);
                atom_fg_strs[it->second] += std::to_string(idx) + "|";
            }
        }
    }

    // For each atom, store the list of functional group indexes in property
    // "FunctionalGroups"
    RDKit::ROMol::AtomIterator ai;
    for (ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai) {
        // Add an additional function group to indicate 'No Functional Groups'
        if (atom_fg_idxs[(*ai)->getIdx()].empty())
            atom_fg_idxs[(*ai)->getIdx()].push_back(num_grps);
        (*ai)->setProp(prop_name, atom_fg_idxs[(*ai)->getIdx()]);
        //(*ai)->setProp("Fg_Idx_Strs", atom_fg_strs[(*ai)->getIdx()]);
    }

    // For each bond figure out if it is between two function groups
    // "FunctionalGroups"
    /*for (auto bi = rwmol->beginBonds(); bi != rwmol->endBonds(); ++bi) {
        std::string f0_fg_idx_str;
        std::string f1_fg_idx_str;
        (*bi)->getBeginAtom()->getProp("Fg_Idx_Strs", f0_fg_idx_str);
        (*bi)->getEndAtom()->getProp("Fg_Idx_Strs", f1_fg_idx_str);

        if(f0_fg_idx_str != f1_fg_idx_str)
            (*bi)->setProp("BetweenFG", 1);
        else
            (*bi)->setProp("BetweenFG", 0);
    }*/
}

void FeatureHelper::labelMMFFAtomTypes(RDKit::RWMol *rwmol) {

    // H-H causes exception...so assign to H-C atom type 5 (which is not used
    // anyway)
    if (rwmol->getAtomWithIdx(0)->getSymbol() == "H") {
        rwmol->getAtomWithIdx(0)->setProp("MMFFAtomType", (int) 5);
        return;
    }
    // For each atom, will store the result in prop "MMFFAtomType"
    RDKit::MMFF::MMFFMolProperties molprop(*rwmol);

    RDKit::ROMol::AtomIterator ai;
    for (ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai) {
        uint8_t atomtype = molprop.getMMFFAtomType((*ai)->getIdx());
        (*ai)->setProp("MMFFAtomType", (int) atomtype);
    }
}

void FeatureHelper::labelAromatics(RDKit::RWMol *rwmol) {

    // Set bond aromaticity information
    RDKit::ROMol::BondIterator bi;
    for (bi = rwmol->beginBonds(); bi != rwmol->endBonds(); ++bi) {
        int aromatic = (*bi)->getIsAromatic();
        (*bi)->setProp("InAromaticRing", aromatic);
        (*bi)->setProp("InDblAromaticRing", 0);
    }

    // Check for any double-aromatic systems
    RDKit::MolOps::findSSSR(*rwmol);
    RDKit::RingInfo *rinfo = rwmol->getRingInfo();
    std::vector<int> double_aromatic_idxs;
    for (unsigned int i = 0; i < rwmol->getNumBonds(); i++) {
        if (rinfo->numBondRings(i) <= 1)
            continue;
        RDKit::Bond *bond = rwmol->getBondWithIdx(i);
        if (bond->getIsAromatic())
            double_aromatic_idxs.push_back(i);
    }

    // If any are found, label all the bonds within them
    if (double_aromatic_idxs.empty())
        return;

    // Consider each ring...
    RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
    auto bit = brings.begin();
    for (; bit != brings.end(); ++bit) {

        // Check for a double aromatic bond within the ring
        bool hasDblArom = false;
        RDKit::RingInfo::INT_VECT::iterator it;
        for (it = bit->begin(); it != bit->end(); ++it) {
            auto ii = double_aromatic_idxs.begin();
            for (; ii != double_aromatic_idxs.end(); ++ii)
                if (*ii == *it)
                    hasDblArom = true;
            if (hasDblArom)
                break;
        }

        // If one exists, label all bonds in the ring
        if (!hasDblArom)
            continue;
        for (it = bit->begin(); it != bit->end(); ++it) {
            RDKit::Bond *bond = rwmol->getBondWithIdx(*it);
            bond->setProp("InDblAromaticRing", 1);
        }
    }
}

void FeatureHelper::labelOriginalMasses(RDKit::RWMol *rwmol) {
    RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
    RDKit::ROMol::AtomIterator ai;
    for (ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai) {
        double mass = 0.0;
        std::string symbol = (*ai)->getSymbol();
        mass += pt->getMostCommonIsotopeMass(symbol);
        mass += (*ai)->getTotalNumHs() * pt->getMostCommonIsotopeMass("H");
        (*ai)->setProp("OriginalMass", mass);
    }
}

void FeatureHelper::labelAtomsWithLonePairs(RDKit::RWMol *rwmol) {
    RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
    RDKit::ROMol::AtomIterator ai;
    RDKit::MolOps::findSSSR(*rwmol);
    
    for (ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai) {
        std::string symbol = (*ai)->getSymbol();
        int nouter = pt->getNouterElecs(symbol.c_str());
        int def_val = pt->getDefaultValence(symbol.c_str());
        (*ai)->setProp(
                "HasLP",
                (int) (nouter > def_val && def_val != 1 &&
                       def_val != -1)); // Allow O,N,S,P..but not C, Halogens, Metals,.
    }
}

int FeatureHelper::getBondTypeAsInt(RDKit::Bond *bond) {
    if (bond->getIsAromatic()) {
        return 4;
    } else if (bond->getIsConjugated()) {
        return 5;
    } else {
        return (int) (bond->getBondTypeAsDouble());
    }
}

void FeatureHelper::labelOriginalBondTypes(RDKit::RWMol *rwmol) {
    RDKit::ROMol::BondIterator bi;
    for (bi = rwmol->beginBonds(); bi != rwmol->endBonds(); ++bi) {
        //(*bi)->setProp("OrigBondTypeRaw", (int)((*bi)->getBondType()));
        int bondType = getBondTypeAsInt(*bi);
        (*bi)->setProp("OrigBondType", bondType);
    }
}