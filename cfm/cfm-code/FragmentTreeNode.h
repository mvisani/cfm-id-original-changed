/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentTreeNode.h
#
# Description: 	Contains Break and FragmentTreeNode classes, for use in
#				the fragment generation process.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FRAG_TREE_NODE_H__
#define __FRAG_TREE_NODE_H__

#include "Util.h"
#include "Feature.h"
#include "Features/FeatureHelper.h"

// Class for storing information about a particular break
class Break {
public:
    // Constructor for a Hydrogen only break
    Break()
            : h_only_break(true), ring_break(false), ionic_break(false), bond_idx(-1),
              ring_idx(-1), num_ionic_frag_allocations(1),
              ionic_idx(-1), is_cycliaztion(false) {};

    // Constructor for Break (ionic or standard or ring)
    Break(int a_bond_or_ionic_atom_idx, bool a_ionic_break,
             bool is_ring_break, int a_ring_idx, int a_num_ionic_frag_allocations, bool is_cycliaztion)
            : h_only_break(false), ring_break(is_ring_break), ionic_break(a_ionic_break),
              bond_idx(-1), ring_idx(-1), ionic_idx(-1),
              num_ionic_frag_allocations(a_num_ionic_frag_allocations), is_cycliaztion(is_cycliaztion) {

        if (is_ring_break)
            ring_idx = a_ring_idx;
        if (ionic_break)
            ionic_idx = a_bond_or_ionic_atom_idx;
        else
            bond_idx = a_bond_or_ionic_atom_idx;
    };

    // Access functions (note: bools convert to ints for use in RDKit properties)
    int getBondIdx() const { return bond_idx; };

    int isRingBreak() const { return ring_break; };

    int isIonicBreak() const { return ionic_break; };

    int isHydrogenOnlyBreak() const { return h_only_break; };

    int getIonicIdx() const { return ionic_idx; };

    int getRingIdx() const { return ring_idx; };

    int getNumIonicFragAllocations() const { return num_ionic_frag_allocations; };

    bool isCycliaztion() const { return is_cycliaztion; };

private:
    int bond_idx;    // Indexes of the broken bond(s)

    bool ring_break; // Flag indicating ring break
    int ring_idx; // Index of the ring that is broken

    bool ionic_break;
    int ionic_idx;
    int num_ionic_frag_allocations;

    bool h_only_break;
    bool is_cycliaztion;
};

// Class for generating fragments via the systematic bond disconnection approach
// with extra checks for hydrogen and bond allocations via MILP.
class FragmentTreeNode {
public:
    // Constructors
    FragmentTreeNode(FeatureHelper *a_fh) : fh(a_fh) {};

    FragmentTreeNode(romol_ptr_t an_ion, romol_ptr_t a_nl, int a_ion_free_ep, int a_depth, FeatureHelper *a_fh,
                         std::vector<int> &a_e_loc, bool is_intermediate = false)
            : ion(an_ion), nl(a_nl), ion_free_epairs(a_ion_free_ep), depth(a_depth),
              fh(a_fh), e_loc(a_e_loc), is_intermediate(is_intermediate) {
        labelIonProperties();
    };

    // used for create start node
    FragmentTreeNode(romol_ptr_t an_ion, int a_ion_free_ep, int a_depth,
                     FeatureHelper *a_fh, std::vector<int> &a_e_loc)
            : ion(an_ion), ion_free_epairs(a_ion_free_ep), depth(a_depth), fh(a_fh),
              e_loc(a_e_loc), is_intermediate(false)  {
        labelIonProperties();
    };

    std::vector<FragmentTreeNode> children;
    romol_ptr_t ion; // The ion fragment
    romol_ptr_t nl; // The neutral loss resulting in this ion
    int ion_free_epairs; // The number of free electron pairs in the ion
    int depth; // The depth of the tree at which the node was

    // Generate all possible breaks of the ion
    // num_rbreak_nrbonds is used to reduce mount of ring breaks
    // if a ion has less than give none ring break bonds, a ring break can happen
    // not the best idea ever
    void generateBreaks(std::vector<Break> &breaks, bool include_H_only_loss, bool include_cyclization);

    // Record a break in the properties of the ion
    void applyBreak(Break &brk, int ionic_allocation_idx);

    // Undo any changes made during applyBreak
    void undoBreak(Break &brk, int ionic_allocation_idx);

    // For an already applied break, generate the possible child fragments
    // and add them to the children field in the node
    void generateChildrenOfBreak(Break &brk);

    void setTmpTheta(double val, int energy) {
        if (tmp_thetas.size() < energy + 1)
            tmp_thetas.resize(energy + 1);
        tmp_thetas[energy] = val;
    };

    double getTmpTheta(int energy) const { return tmp_thetas[energy]; };

    const std::vector<double> *getAllTmpThetas() const { return &tmp_thetas; };

    // Static Utility functions:

    // For a given molecule, find a good place to allocate the charge on the
    // specified fragment
    //: returns atom index of the charge, and atom index of the radical (or -1 if
    //either are unallocated)
    static std::pair<int, int> findChargeLocation(RDKit::RWMol &rwmol,
                                                  int charge_side,
                                                  int radical_side,
                                                  bool is_negative);

    static boost::tuple<bool, bool, bool>
    findAlreadyChargedOrSplitCharge(boost::tuple<int, int, int> &pidx_nidx_ridx,
                                    RDKit::RWMol &rwmol, int charge_side,
                                    int radical_side, bool is_negative);

    // Apply the charge and radical (if specified) at the provided atom index, and
    // also remove them
    static void assignChargeAndRadical(RDKit::RWMol &rwmol, int charge_idx,
                                       int radical_idx, bool is_negative);

    static void undoChargeAndRadical(RDKit::RWMol &rwmol, int charge_idx,
                                     int radical_idx, bool is_negative);

    static void
    undoAlreadyChargedOrSplitCharge(RDKit::RWMol &rwmol,
                                    boost::tuple<int, int, int> &pidx_nidx_ridx);

    bool isIntermediate() const { return is_intermediate; };
    
    bool isCyclization() const { return is_cyclization; };

    void setCyclization(bool flag) { this->is_cyclization = flag; };

private:
    // Helper class for feature labels that need to be added during fragment graph
    // computation
    FeatureHelper *fh;

    // Storage for where the original electrons were
    std::vector<int> e_loc;

    // Temporary storage for the current theta value in each node
    std::vector<double> tmp_thetas;

    // Helper functions:
    // Used to produce FragIdx labels for a broken molecule
    // assume all FragIdx set to 0 already
    // set atom's FragIdx to 1 recursively until we meet the broken bond
    void allocatedCtdToFragment(RDKit::ROMol *romol, RDKit::Atom *atom);

    // Creates all the details and adds the children (with charge on either side)
    // of the node. -- e_f0 specifies the number of epairs to assigne to F0,
    // e_to_allocate is the total number to allocate, and output_bmax specifies
    // where the epairs can go.
    int addChild(int e_f0, int e_to_allocate, std::vector<int> &output_bmax,
                  Break &brk, int charge_frag);

    // Utility function for counting the number of electron pairs allocated to the
    // original fragments
    std::pair<int, int> computeOrigFreeElectronsPerFrag();

    // Utility function to store various properties of the break in the generated
    // neutral loss
    void labelBreakPropertiesInNL(romol_ptr_t &current_nl, romol_ptr_t &parent_ion,
                                  Break &brk);

    // Utility function to store extra  properties of the break in the generated iom
    void labelExtraBreakPropertiesInIon(romol_ptr_t &current_ion, romol_ptr_t &parent_ion,
                                  Break &brk);


    // Utility function to label properties of the ion: is it negative? is it a
    // radical?
    void labelIonProperties();

    void createChildIonElectronLocRecord(std::vector<int> &child_e_loc,
                                         romol_ptr_t childmol);

    // Static Utility functions:
    static int findAtomChargeLocationNSOC(RDKit::RWMol &rwmol, int charge_side,
                                          bool is_negative);

    static int countNumIonicFragments(
            const RDKit::ROMol
            *romol); // Requires IonicFragmentCharge labels to be present
    static int computeNumIonicAlloc(int num_ionic_fragments);

    static void recordOrigAtomIdxs(RDKit::RWMol &rwmol);

    // Flags
    bool is_intermediate = false;
    bool is_cyclization = false;
};

#endif // __FRAG_TREE_NODE_H__