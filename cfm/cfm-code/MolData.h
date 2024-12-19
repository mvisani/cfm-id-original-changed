/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# MolData.h
#
# Description: 	Class to hold the input data belonging to a molecule:
#					 - An ID and smiles/inchi
#					 - (optional) A cross-validation group
identifier
#					 - (optional) A computed fragmentation
graph
#					 - (optional) A computed set of
features corresponding to that graph
#					 - (optional) A computed set of theta
values for that graph
#				     - (optional) A computed set of transition
probabilities using those thetas.
#					 - (optional) A set of spectra
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __MOLDATA_H__
#define __MOLDATA_H__

#include <DataStructs/ExplicitBitVect.h>

#include "Config.h"
#include "Feature.h"
#include "FragmentGraph.h"
#include "FragmentGraphGenerator.h"
#include "Message.h"
#include "MspReader.h"
#include "NNParam.h"
#include "Param.h"
#include "Spectrum.h"

class MolData {
public:
	MolData(std::string &an_id, std::string &an_smiles_or_inchi, int a_group, config_t *a_cfg)
	    : id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(a_group), graph_computed(false),
	      ev_graph_computed(false), cfg(a_cfg) {};

	MolData(const char *an_id, const char *an_smiles_or_inchi, config_t *a_cfg)
	    : id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(0), graph_computed(false), ev_graph_computed(false),
	      cfg(a_cfg) {};

	// Access functions
	const FragmentGraph *getFragmentGraph() const { return fg; };

	bool hasComputedGraph() const { return graph_computed; };

	const EvidenceFragmentGraph *getEvidenceFragmentGraph() const { return ev_fg; };

	const Spectrum *getSpectrum(int energy) const { return &(spectra[energy]); };

	const Spectrum *getOrigSpectrum(int energy) const { return &(orig_spectra[energy]); };

	const std::vector<Spectrum> *getSpectra() const { return &spectra; };

	const Spectrum *getPredictedSpectrum(int energy) const { return &(predicted_spectra[energy]); };

	unsigned int getNumSpectra() const { return spectra.size(); };

	unsigned int getNumPredictedSpectra() const { return predicted_spectra.size(); };

	const FeatureVector *getFeatureVectorForIdx(int index) const {
		return getTransitionAtIdx(index)->getFeatureVector();
	};

	double getThetaForIdx(int energy, int index) const { return thetas[energy][index]; };

	double getLogTransitionProbForIdx(int energy, int index) const { return log_probs[energy][index]; };

	double getLogPersistenceProbForIdx(int energy, int index) const {
		return log_probs[energy][fg->getNumTransitions() + index];
	};

	int getGroup() const { return group; };

	void setGroup(int val) { group = val; };

	std::string getId() const { return id; };

	std::string getSmilesOrInchi() const { return smiles_or_inchi; };

	double getMolecularWeight() const;

	double getParentIonMass() const;

	int getIonizationMode() const { return cfg->ionization_mode; };

	void readInSpectraFromFile(const std::string &filename, bool readToPredicted = false);

	void readInSpectraFromMSP(MspReader &msp, bool readToPredicted = false);

	void cleanSpectra(double abs_tol, double ppm_tol);

	void freeSpectra() {
		std::vector<Spectrum>().swap(spectra);
		std::vector<Spectrum>().swap(predicted_spectra);
	};

	// Spectrum Related Functions
	std::string removePeaksWithNoFragment(double abs_tol, double ppm_tol);

	bool hasEmptySpectrum(int energy_level = -1) const;

	void writePredictedSpectraToFile(std::string &filename);

	void writePredictedSpectraToMspFileStream(std::ostream &out);

	void writePredictedSpectraToMgfFileStream(std::ostream &out);

	void writeFullEnumerationSpectrumToFile(std::string &filename);

	void writeFullEnumerationSpectrumToMspFileStream(std::ostream &out);

	void outputSpectra(std::ostream &out, const char *spec_type, bool do_annotate = false, bool add_version = true);

	// More memory efficient alternative to calling computeFragmentGraph and
	// then computeFeatureVectors with deleteMols = true
	void computeFragmentGraphAndReplaceMolsWithFVs(FeatureCalculator *fc, bool retain_smiles = false);

	// Save/load state functions
	void readInFVFragmentGraph(std::string &fv_filename);

	void readInFVFragmentGraphFromStream(std::istream &ifs);

	void writeFVFragmentGraph(std::string &fv_filename);

	void writeFVFragmentGraphToStream(std::ofstream &out);

	// Replaces computeFragmentGraph, computeFeatureVectors and
	// computeTransitionThetas  below (delteMols = true), pruning according to
	// prob_thresh_for_prune value.
	void computeLikelyFragmentGraphAndSetThetas(LikelyFragmentGraphGenerator &fgen, bool retain_smiles);

	// Note that the following should be called in this order
	// since each one assumes all previous have already been called.E
	void computeFragmentGraph(FeatureCalculator *fc);

	void computeTransitionThetas(Param &param);

	void computeLogTransitionProbabilities();

	// compute predicted Spectra
	// if engry < -1 , compute all  Spectra
	void computePredictedSpectra(Param &param, bool use_existing_thetas, int energy_level, int min_peaks, int max_peaks,
	                             double perc_thresh, double min_relative_intensity, int quantise_peaks_decimal_place,
	                             bool log_to_linear);

	void postprocessPredictedSpectra(double perc_thresh, int min_peaks, int max_peaks,
	                                 double min_relative_intensity_prec, int quantise_peaks_decimal_place);

	void quantisePredictedSpectra(int num_dec_places);

	void quantiseMeasuredSpectra(int num_dec_places);

	// Function to compute a much reduced fragment graph containing only those
	// fragmentations as actually occur in the spectra, based on a computed set of
	// beliefs  thresholding inclusion in the graph by the provided belief_thresh
	// value (log domain)
	void computeEvidenceFragmentGraph(beliefs_t *beliefs, double log_belief_thresh);

	void annotatePeaks(double abs_tol, double ppm_tol, bool prune_deadends = true);

	void getSampledTransitionIdsWeightedRandomWalk(std::set<int> &selected_ids, int max_num_iter, int energy,
	                                               double explore_weight);

	void getSampledTransitionIdsRandomWalk(std::set<int> &selected_ids, int max_selection);

	void getSampledTransitionIdUsingDiffMapBFS(std::set<int> &selected_ids, std::set<unsigned int> &selected_weights);

	void getSampledTransitionIdUsingDiffMapCA(std::set<int> &selected_ids, std::set<unsigned int> &selected_weights);

	void getRandomSampledTransitions(std::set<int> &selected_ids, int max_selection);

	unsigned int getNumTransitions() const { return fg->getNumTransitions(); };

	unsigned int getNumFragments() const { return fg->getNumFragments(); };

	const TransitionPtr getTransitionAtIdx(int index) const { return fg->getTransitionAtIdx(index); };

	virtual const Fragment *getFragmentAtIdx(int index) const { return fg->getFragmentAtIdx(index); };

	const tmap_t *getFromIdTMap() const { return fg->getFromIdTMap(); };

	const tmap_t *getToIdTMap() const { return fg->getToIdTMap(); };

	void writeFullGraph(std::ostream &out) const { fg->writeFullGraph(out); };

	void writeFragmentsOnly(std::ostream &out) const { fg->writeFragmentsOnly(out); }

	void writeFragmentsOnlyForIds(std::ostream &out, std::set<int> &ids) const {
		fg->writeFragmentsOnlyForIds(out, ids);
	}

	bool hasIsotopesIncluded() const { return fg->hasIsotopesIncluded(); }

	int getFGHeight() const { return fg->getHeight(); };

	// Return a listed of select weights in Fixed Point INT
	// int_weight = std::round(weight * 1e-5)
	void getSelectedMasses(std::set<unsigned int> &selected_weights, int energry_level);

	void computeMergedPrediction();

	const Spectrum *getMergedPrediction() { return m_merged_predicted_spectra; };

	void convertSpectraToLogScale();

	void convertSpectraToLinearScale();

	~MolData();

protected: // These items are protected rather than private for access during tests.
	int group;
	std::string id;
	std::string smiles_or_inchi;
	FragmentGraph *fg            = nullptr;
	EvidenceFragmentGraph *ev_fg = nullptr;
	bool graph_computed;
	bool ev_graph_computed;
	// spectra , which will be pruned during training
	std::vector<Spectrum> spectra;
	// orig copy of spectra
	std::vector<Spectrum> orig_spectra;
	// predicted spectra
	std::vector<Spectrum> predicted_spectra;

	// merged predicted spectra
	// used in casmi
	Spectrum *m_merged_predicted_spectra = nullptr;

	// std::vector<FeatureVector *> fvs;
	std::vector<std::vector<double>> thetas;
	std::vector<std::vector<double>> log_probs;
	config_t *cfg = nullptr;

	// General utilty functions
	void computeGraphWithGenerator(FragmentGraphGenerator &fgen);

	void getEnumerationSpectraMasses(std::vector<double> &output_masses);

	void translatePeaksFromMsgToSpectra(Spectrum &out_spec, Message *msg);

	void translatePeaksFromMsgToSpectraWithIsotopes(Spectrum &out_spec, Message *msg);

	void computeFragmentEvidenceValues(std::vector<double> &evidence, int frag_idx, const beliefs_t *beliefs);

	void createSpeactraSingleEnergry(unsigned int energy_level);
};

#endif // __MOLDATA_H__
