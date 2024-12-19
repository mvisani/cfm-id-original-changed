/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraph.h
#
# Description: FragmentGraph class for holding the results of a generated
#               fragment graph.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FRAGTREE_H__
#define __FRAGTREE_H__

#include <GraphMol/ROMol.h>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <vector>

#include "Feature.h"
#include "Features/FeatureHelper.h"
#include "FragmentTreeNode.h"
#include "Isotope.h"
#include "Util.h"

typedef std::vector<std::vector<int>> tmap_t;
// Class for storing the base fragment state for our model
class Fragment {

public:
	// Constructor, store the ion smiles and a reduced smiles since the ion is
	// not needed and takes more space.
	Fragment() {};

	Fragment(std::string &a_ion_smiles, std::string &a_reduced_smiles, int an_id, double a_mass, bool is_intermediate,
	         bool is_cyclization)
	    : id(an_id), ion_smiles(a_ion_smiles), reduced_smiles(a_reduced_smiles), mass(a_mass), depth(-1),
	      is_intermediate(is_intermediate), is_cyclization(is_cyclization) {};

	Fragment(std::string &a_ion_smiles, std::string &a_reduced_smiles, int an_id, double a_mass,
	         Spectrum &a_isotope_spec, bool is_intermediate, bool is_cyclization)
	    : id(an_id), ion_smiles(a_ion_smiles), reduced_smiles(a_reduced_smiles), mass(a_mass),
	      isotope_spectrum(a_isotope_spec), depth(-1), is_intermediate(is_intermediate),
	      is_cyclization(is_cyclization) {};

	Fragment(const Fragment &a_fragment, int an_id)
	    : id(an_id), ion_smiles(*a_fragment.getIonSmiles()), reduced_smiles(*a_fragment.getReducedSmiles()),
	      mass(a_fragment.getMass()), depth(-1), is_intermediate(a_fragment.isIntermediate()),
	      is_cyclization(a_fragment.isCyclization()), isotope_spectrum(*a_fragment.getIsotopeSpectrum()) {};

	// Access Functions
	double getMass() const { return mass; };

	const Spectrum *getIsotopeSpectrum() const { return &isotope_spectrum; };

	int getId() const { return id; };

	void setId(int an_id) { id = an_id; };

	void setIntermediate(bool flag) { is_intermediate = flag; };

	bool isIntermediate() const { return is_intermediate; };

	void setCyclization(bool flag) { is_cyclization = flag; };

	bool isCyclization() const { return is_cyclization; };

	const std::string *getReducedSmiles() const { return &reduced_smiles; };

	const std::string *getIonSmiles() const { return &ion_smiles; };

	void clearSmiles() {
		reduced_smiles = std::string();
		ion_smiles     = std::string();
	};

	void setDepth(int a_depth) { depth = a_depth; };

	int getDepth() const { return depth; };

protected:
	int id;
	// RootedROMol *ion = nullptr;
	std::string reduced_smiles; // Reduced version of the smiles string (just backbone)
	std::string ion_smiles;     // Full ion smiles (for writing out if called for)
	double mass;
	Spectrum isotope_spectrum;
	int depth; // Depth -1 means hasn't been set yet.
	bool is_intermediate = false;
	bool is_cyclization  = false;
};

class EvidenceFragment : public Fragment {
public:
	EvidenceFragment(const Fragment &a_fragment, int an_id, const std::vector<double> &a_evidence)
	    : Fragment(a_fragment, an_id), evidence(a_evidence) {};

	EvidenceFragment(const EvidenceFragment &a_fragment, int an_id)
	    : Fragment(a_fragment, an_id), evidence(a_fragment.evidence) {};

	double getEvidence(int energy) const { return evidence[energy]; };

protected:
	std::vector<double> evidence; // For recording the belief that this fragment
	                              // occurs at each energy level (in an evidence
	                              // fragment graph - used in peak annotation)
};

// Class for storing possible transitions between fragments
class Transition {

public:
	// Default constructor
	Transition() {};

	// Basic constructor
	Transition(int a_from_id, int a_to_id, const RootedROMol &a_nl, const RootedROMol &an_ion);

	// Alternative constructor that finds the root atoms and sets the
	// root pointers appropriately
	Transition(int a_from_id, int a_to_id, const romol_ptr_t &a_nl, const romol_ptr_t &an_ion);

	// Direct constructor that bipasses the mols altogether and directly sets the
	// nl_smiles
	Transition(int a_from_id, int a_to_id, const std::string *a_nl_smiles)
	    : from_id(a_from_id), to_id(a_to_id), nl_smiles(*a_nl_smiles) {};

	~Transition() {
		if (feature_vector != nullptr) delete feature_vector;
	};

	// Access Functions
	int getFromId() const { return from_id; };

	int getToId() const { return to_id; };

	void setFromId(const int id) { from_id = id; };

	void setToId(const int id) { to_id = id; };

	const std::string *getNLSmiles() const { return &nl_smiles; };

	const RootedROMol *getNeutralLoss() const { return &nl; };

	void deleteNeutralLoss() {
		nl.mol.reset();
		nl = RootedROMol();
	};

	const RootedROMol *getIon() const { return &ion; };

	void deleteIon() {
		ion.mol.reset();
		ion = RootedROMol();
	};

	FeatureVector *getFeatureVector() const { return feature_vector; };

	void setFeatureVector(FeatureVector *an_fv_ptr) { feature_vector = an_fv_ptr; };

	const std::vector<double> *getTmpThetas() const { return &tmp_thetas; };

	void setTmpThetas(const std::vector<double> *a_thetas) { tmp_thetas = *a_thetas; };

	bool isDuplicate() { return is_duplicate; }

	void createdDuplication(const Transition &old) {
		from_id   = old.from_id;
		to_id     = old.to_id;
		nl_smiles = old.nl_smiles;
		if (old.feature_vector != nullptr) feature_vector = new FeatureVector(*old.feature_vector);
		is_duplicate = true;
		tmp_thetas   = old.tmp_thetas;
	};

private:
	int from_id;
	int to_id;
	// we should ONLY keep nl smiles
	// since ion could be re-used from an older fragments
	std::string nl_smiles;
	RootedROMol nl;

	// This is a BUG  .... need fix
	// TODO:fix this
	// if we are using ion
	// this should be a re-used ion
	RootedROMol ion;

	// We store the ion on the transition to
	// allow for different roots - the fragment stores
	// only an unrooted shared pointer.
	FeatureVector *feature_vector = nullptr; // storage for the feature vector pointer

	// while we compute the fragment graph (don't use this
	// directly - it will be moved up into the MolData)
	std::vector<double> tmp_thetas;
	// Temporary storage for the theta values
	// while we compute the likely fragment graph
	//(as above, don't use directly, will be moved)

	// a flag
	bool is_duplicate = false;
};

typedef std::shared_ptr<Transition> TransitionPtr;

class FragmentGraph {
public:
	FragmentGraph()
	    : include_isotopes(false), allow_frag_detours(true), include_h_losses(true),
	      include_h_losses_precursor_only(false), allow_cyclization(false) {};

	FragmentGraph(config_t *cfg)
	    : include_isotopes(cfg->include_isotopes), allow_frag_detours(cfg->allow_frag_detours),
	      include_h_losses(cfg->include_h_losses),
	      include_h_losses_precursor_only(cfg->include_precursor_h_losses_only),
	      allow_cyclization(cfg->allow_cyclization) {
		if (include_isotopes) isotope = new IsotopeCalculator(cfg->isotope_thresh, cfg->isotope_pattern_file);
	};

	~FragmentGraph() {
		if (include_isotopes) delete isotope;
		for (auto &fragment : fragments) { delete fragment; }
	};

	// Add a fragment node to the graph (should be the only way to modify the
	// graph)
	//	-- Add a fragment, or return an id, if it already exists
	//  -- Add a transition to this fragment, based on the provided parent
	//  fragment id
	//  -- Update the relevant tmaps
	// Note: If parentid < 0, assumes starting ion, so adds extra H to mass, and
	//       doesn't add transition.
	int addToGraph(const FragmentTreeNode &node, int parentid);

	// As for previous function, but delete the mols in the transition and compute
	// and store a feature vector instead
	int addToGraphAndReplaceMolWithFV(const FragmentTreeNode &node, int parent_frag_id, FeatureCalculator *fc);

	// As for previous function, but don't store the mols in the transition and
	// insert the pre-computed thetas instead
	int addToGraphWithThetas(const FragmentTreeNode &node, const std::vector<double> *thetas, int parent_frag_id);

	// Write the Fragments only to file (formerly the backtrack output - without
	// extra details)
	virtual void writeFragmentsOnly(std::ostream &out) const;

	virtual void writeFragmentsOnlyForIds(std::ostream &out, std::set<int> &ids) const;

	// Write the FragmentGraph to file (formerly the transition output - without
	// feature details)
	void writeFullGraph(std::ostream &out) const;

	// Write the fragment graph - no smiles, just ids, fragment masses and feature
	// vectors for the transitions.
	void writeFeatureVectorGraph(std::ostream &out, bool include_isotopes) const;

	void readFeatureVectorGraph(std::istream &out);
	;
	;

	void addFeatureVectorAtIdx(int index, FeatureVector *feature_vector) const {
		transitions[index]->setFeatureVector(feature_vector);
	};

	bool hasIsotopesIncluded() const { return include_isotopes; };

	double getIsotopeThresh() const { return isotope->getIntensityThresh(); };

	bool includesHLosses() const { return include_h_losses; };

	bool includesHLossesPrecursorOnly() const { return include_h_losses_precursor_only; };

	bool allowCyclization() const { return allow_cyclization; };

	void clearAllSmiles();

	// For current graph in use
	virtual unsigned int getNumTransitions() const { return transitions.size(); };

	virtual unsigned int getNumFragments() const { return fragments.size(); };

	virtual const TransitionPtr getTransitionAtIdx(int index) const { return transitions[index]; };

	const Fragment *getFragmentAtIdx(int index) const { return fragments[index]; };

	const tmap_t *getFromIdTMap() const { return &from_id_tmap; };

	const tmap_t *getToIdTMap() const { return &to_id_tmap; };

	// Function to remove detour transitions from the graph (used if
	// !cfg.allow_frag_detours)
	void removeDetours();

	// Get a list of transitions ids , with weighted prob
	// Function do some not so random selection
	void getSampledTransitionIdsWeightedRandomWalk(std::set<int> &selected_ids, int max_num_iter,
	                                               std::vector<double> &thetas, double explore_weight);

	// Get a list of transitions ids , with random walk
	void getSampledTransitionIdsRandomWalk(std::set<int> &selected_ids, int max_selection);

	// Get a list of transitions ids , with random walk
	void getSampledTransitionIdsDiffMapChildOnly(std::set<int> &selected_ids, std::set<unsigned int> &selected_weights);

	void getSampledTransitionIdsDiffMap(std::set<int> &selected_ids, std::set<unsigned int> &selected_weights);

	void getRandomSampledTransitions(std::set<int> &selected_trans_id, int max_selection);

	int getHeight() { return depth; };

	int setHeight(int d) { return depth = d; };

protected:
	std::vector<Fragment *> fragments;
	std::vector<TransitionPtr> transitions;
	int depth = 0;

	tmap_t from_id_tmap; // Mapping between from_id and transitions with that from_id
	tmap_t to_id_tmap;   // Mapping between to_id and transitions with that to_id

	void getSampledTransitionIdsWeightDiffChildOnly(std::set<unsigned int> &selected_weights, std::set<int> &visited,
	                                                int frag_id, std::vector<int> &path, std::set<int> &selected_ids);

	void getSampledTransitionIdsWeightDiffs(std::set<unsigned int> &selected_weights, std::set<int> &visited,
	                                        int frag_id, std::vector<std::pair<int, int>> &path,
	                                        std::map<int, std::vector<int>> &trans_to_interest_frags_map);

	bool is_match(std::set<unsigned int> &weights, double mass) const;

	bool include_isotopes;
	IsotopeCalculator *isotope;
	bool allow_frag_detours;
	bool include_h_losses;
	bool include_h_losses_precursor_only;
	bool allow_cyclization;

	// Mapping from rounded mass to list of fragment ids,
	// to enable fast check for existing fragments
	std::map<double, std::vector<int>> frag_mass_lookup;

	// Find the id for an existing fragment that matches the input ion and mass
	// or create a new fragment in the case where no such fragment is found
	int addFragmentOrFetchExistingId(romol_ptr_t ion, double mass, bool is_intermediate, bool is_cyclization);

	// Determine if the two fragments match - assumes the masses have already
	// been checked to be roughly the same, now check reduced structure.
	bool areMatching(RDKit::ROMol *f1_reduced_ion, RDKit::ROMol *f2_reduced_ion);

	// Make all bonds single and remove any charge - used to compare structure
	// alone
	void reduceMol(RDKit::RWMol &rwmol);

	// Find the id for an existing transition that matches the input ids
	// or -1 in the case where no such transition is found
	int findMatchingTransition(int from_id, int to_id);
};

class EvidenceFragmentGraph : public FragmentGraph {
public:
	EvidenceFragmentGraph(config_t *cfg) : FragmentGraph(cfg) {};

	const EvidenceFragment *getFragmentAtIdx(int index) const { return &(fragments[index]); };

	void writeFragmentsOnly(std::ostream &out) const override;

	void writeFullGraph(std::ostream &out) const;

	void writeFragmentsOnlyForIds(std::ostream &out, std::set<int> &ids) const override;

	// Alternative graph builder used to build fragments from fragments and
	// transitions of an existing graph
	// Note: no checks for duplicates, assumes the added fragments have already
	// been filtered.
	int addToGraphDirectNoCheck(const EvidenceFragment &fragment, const Transition *transition, int parentid);

	void addTransition(int from_id, int to_id, const std::string *nl_smiles);

	// Utility functions used in annotation
	bool fragmentIsRedundant(unsigned int fidx, std::vector<int> &annotated_flags, std::vector<int> &direct_flags,
	                         std::vector<bool> &visited) const;

	void setFlagsForDirectPaths(std::vector<int> &direct_flags, unsigned int fidx,
	                            std::vector<int> &annotated_flags) const;

	// Access functions, override base functions
	unsigned int getNumTransitions() const override { return transitions.size(); };

	unsigned int getNumFragments() const override { return fragments.size(); };

	const TransitionPtr getTransitionAtIdx(int index) const override { return transitions[index]; };

private:
	std::vector<EvidenceFragment> fragments;
};

#endif // __FRAGTREE_H__
