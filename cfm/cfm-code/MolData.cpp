/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# MolData.cpp
#
# Description: 	Class to hold the input data belonging to a molecule:
#					 - An ID and smiles/inchi
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

#include "MolData.h"
#include "Comparators.h"
#include "Inference.h"
#include "Version.h"

#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/inchi.h>
#include <boost/algorithm/string.hpp>
#include <regex>
#include <string>

double MolData::getMolecularWeight() const {
	romol_ptr_t mol = createMolPtr(smiles_or_inchi.c_str());
	return getMonoIsotopicMass(mol);
}

void MolData::readInFVFragmentGraphFromStream(std::istream &ifs) {

	FragmentGraphGenerator fgen;
	fg = fgen.createNewGraph(cfg);
	fg->readFeatureVectorGraph(ifs);

	graph_computed = true;
}

void MolData::readInFVFragmentGraph(std::string &fv_filename) {

	std::ifstream ifs(fv_filename.c_str(), std::ifstream::in | std::ios::binary);

	if (!ifs) std::cerr << "Could not open file " << fv_filename << std::endl;
	readInFVFragmentGraphFromStream(ifs);
}

void MolData::writeFVFragmentGraphToStream(std::ofstream &out) {
	fg->writeFeatureVectorGraph(out, cfg->include_isotopes);
}

void MolData::writeFVFragmentGraph(std::string &fv_filename) {

	std::ofstream out;
	out.open(fv_filename.c_str(), std::ios::out | std::ios::binary);
	if (!out.is_open()) {
		std::cerr << "Warning: Trouble opening output fv fragment graph file: " << fv_filename << std::endl;
	} else {
		writeFVFragmentGraphToStream(out);
	}
}

void MolData::convertSpectraToLogScale() {
	for (auto &spectrum : spectra) spectrum.convertToLogScale();
}

void MolData::convertSpectraToLinearScale() {
	for (auto &spectrum : spectra) spectrum.convertToLinearScale();
}

void MolData::computeGraphWithGenerator(FragmentGraphGenerator &fgen) {

	try {
		fg                          = fgen.createNewGraph(cfg);
		FragmentTreeNode *startnode = fgen.createStartNode(smiles_or_inchi, cfg->ionization_mode);

		const int root_id = -1;
		fgen.compute(*startnode, cfg->fg_depth, root_id, cfg->max_ring_breaks);

		if (!cfg->allow_frag_detours) fg->removeDetours();
		delete startnode;
		graph_computed = true;
	} catch (std::exception &e) {
		std::cerr << "Error Processing " << smiles_or_inchi << " " << e.what() << std::endl;
		graph_computed = false;
	}
}

void MolData::computeFragmentGraph(FeatureCalculator *fc) {

	FragmentGraphGenerator fgen(fc, false);
	computeGraphWithGenerator(fgen);
}

void MolData::computeFragmentGraphAndReplaceMolsWithFVs(FeatureCalculator *fc, bool retain_smiles) {

	// Compute the fragment graph, replacing the transition molecules with feature
	// vectors
	FragmentGraphGenerator fgen(fc);
	computeGraphWithGenerator(fgen);

	// Delete all the fragment smiles (we only need these while we're computing
	// the graph)
	if (!retain_smiles) fg->clearAllSmiles();
}

void MolData::computeLikelyFragmentGraphAndSetThetas(LikelyFragmentGraphGenerator &fgen, bool retain_smiles) {
	// Compute the fragment graph, replacing the transition molecules with feature
	// vectors
	fg                          = fgen.createNewGraph(cfg);
	FragmentTreeNode *startnode = fgen.createStartNode(smiles_or_inchi, cfg->ionization_mode);
	fgen.compute(*startnode, cfg->fg_depth, -1, 0.0, cfg->max_ring_breaks);
	if (!cfg->allow_frag_detours) fg->removeDetours();

	delete startnode;
	graph_computed = true;

	// Copy all the theta values up into the mol data and delete the tmp thetas
	const unsigned int num_levels = cfg->spectrum_depths.size();
	thetas.resize(num_levels);
	for (unsigned int energy = 0; energy < num_levels; energy++) {
		thetas[energy].resize(fg->getNumTransitions());
		for (unsigned int i = 0; i < fg->getNumTransitions(); i++)
			thetas[energy][i] = (*(fg->getTransitionAtIdx(i)->getTmpThetas()))[energy];
	}

	// Delete all the fragment smiles (we only need these while we're computing
	// the graph)
	if (!retain_smiles) fg->clearAllSmiles();
}

// Function to compute a much reduced fragment graph containing only those
// fragmentations as actually occur in the spectra, based on a computed set of
// beliefs  thresholding inclusion in the graph by the provided belief_thresh
// value (log domain)
void MolData::computeEvidenceFragmentGraph(beliefs_t *beliefs, double log_belief_thresh) {

	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments   = fg->getNumFragments();
	ev_fg                        = new EvidenceFragmentGraph(cfg);

	// Add the root fragment
	Transition null_trans;
	std::vector<int> id_lookup(num_fragments, -1);
	std::vector<double> main_ev;
	computeFragmentEvidenceValues(main_ev, 0, beliefs);
	id_lookup[0] =
	    ev_fg->addToGraphDirectNoCheck(EvidenceFragment(*(fg->getFragmentAtIdx(0)), -1, main_ev), &null_trans, -1);

	// Add the fragments and transitions if the beliefs are high enough
	std::vector<int> t_added_flags(num_transitions, 0);
	for (int depth = 0; depth < beliefs->tn[0].size(); depth++) {

		for (unsigned int i = 0; i < num_transitions; i++) {
			if (t_added_flags[i]) continue;

			if (beliefs->tn[i][depth] >= log_belief_thresh) {
				const TransitionPtr t = fg->getTransitionAtIdx(i);
				const Fragment *f     = fg->getFragmentAtIdx(t->getToId());
				int from_id           = id_lookup[t->getFromId()];
				if (id_lookup[t->getToId()] == -1) {
					std::vector<double> evidence;
					computeFragmentEvidenceValues(evidence, t->getToId(), beliefs);
					int to_id = ev_fg->addToGraphDirectNoCheck(EvidenceFragment(*f, -1, evidence), t.get(), from_id);
					id_lookup[t->getToId()] = to_id;
				} else if (id_lookup[t->getFromId()] != -1)
					ev_fg->addTransition(from_id, id_lookup[t->getToId()], t->getNLSmiles());
				t_added_flags[i] = 1;
			}
		}
	}
	ev_graph_computed = true;
}

void MolData::computeFragmentEvidenceValues(std::vector<double> &evidence, int frag_idx, const beliefs_t *beliefs) {

	const tmap_t *tomap = fg->getToIdTMap();

	evidence.resize(cfg->dv_spectrum_depths.size());
	for (int energy = 0; energy < cfg->dv_spectrum_depths.size(); energy++) {

		// Find out the depth of the spectrum of interest in the beliefs
		int depth = cfg->dv_spectrum_depths[energy] - 1;

		for (int i = 0; i < energy; i++) depth += cfg->dv_spectrum_depths[i];

		// Compute the accumulated belief for the fragment of interest at the depth
		// of interest
		evidence[energy]                    = beliefs->ps[frag_idx][depth];
		std::vector<int>::const_iterator it = (*tomap)[frag_idx].begin();
		for (; it != (*tomap)[frag_idx].end(); ++it)
			evidence[energy] = logAdd(evidence[energy], beliefs->tn[*it][depth]);
	}
}

void MolData::annotatePeaks(double abs_tol, double ppm_tol, bool prune_deadends) {

	if (!ev_graph_computed) {
		std::cout << "Warning: called annotate peaks without first computing "
		             "evidence graph"
		          << std::endl;
		return;
	}

	unsigned int num_fragments = ev_fg->getNumFragments();
	std::vector<int> frag_flags(num_fragments, 0);

	// Assign fragment annotations to peaks
	for (int energy = 0; energy < spectra.size(); energy++) {
		Spectrum::const_iterator it = spectra[energy].begin();
		for (int pkidx = 0; it != spectra[energy].end(); ++it, pkidx++) {

			double mass_tol = getMassTol(abs_tol, ppm_tol, it->mass);
			// int num_fragments = ev_fg->getNumFragments();
			for (int fidx = 0; fidx < num_fragments; fidx++) {

				const EvidenceFragment *f = ev_fg->getFragmentAtIdx(fidx);
				if (ev_fg->hasIsotopesIncluded()) {
					const Spectrum *isotope_spec = f->getIsotopeSpectrum();
					Spectrum::const_iterator itp = isotope_spec->begin();
					for (; itp != isotope_spec->end(); ++itp) {
						if (fabs(itp->mass - it->mass) <= mass_tol) {
							spectra[energy].addPeakAnnotation(pkidx,
							                                  annotation_t(f->getId(), exp(f->getEvidence(energy))));
							frag_flags[fidx] = 1;
						}
					}

				} else if (fabs(f->getMass() - it->mass) <= mass_tol) {
					spectra[energy].addPeakAnnotation(pkidx, annotation_t(f->getId(), exp(f->getEvidence(energy))));
					frag_flags[fidx] = 1;
				}
			}
		}
	}

	// Sort the annotations by evidence score
	for (int energy = 0; energy < spectra.size(); energy++) spectra[energy].sortAndNormalizeAnnotations();

	if (!prune_deadends) return;

	std::vector<int> delete_flags(num_fragments);
	for (unsigned int i = 0; i < num_fragments; i++) {
		// If there is no peak for this fragment, check there is
		// an annotated descendent, else flag for deletion
		std::vector<int> direct_flags(num_fragments);
		std::vector<bool> visited(num_fragments);
		ev_fg->setFlagsForDirectPaths(direct_flags, 0, frag_flags);
		if (ev_fg->fragmentIsRedundant(i, frag_flags, direct_flags, visited))
			delete_flags[i] = 1;
		else
			delete_flags[i] = 0;
	}

	// Rebuild the graph again without deleted fragments
	Transition null_t;
	EvidenceFragmentGraph *new_ev_fg = new EvidenceFragmentGraph(cfg);
	std::vector<int> id_lookup(num_fragments, -1);
	id_lookup[0] = new_ev_fg->addToGraphDirectNoCheck(*ev_fg->getFragmentAtIdx(0), &null_t, -1);
	for (unsigned int i = 0; i < ev_fg->getNumTransitions(); i++) {
		const TransitionPtr t = ev_fg->getTransitionAtIdx(i);
		int fromid            = t->getFromId();
		int toid              = t->getToId();
		if (delete_flags[toid] || delete_flags[fromid]) continue;
		if (id_lookup[toid] == -1)
			id_lookup[toid] =
			    new_ev_fg->addToGraphDirectNoCheck(*ev_fg->getFragmentAtIdx(toid), t.get(), id_lookup[fromid]);
		else
			new_ev_fg->addTransition(id_lookup[fromid], id_lookup[toid], t->getNLSmiles());
	}
	delete ev_fg;
	ev_fg = new_ev_fg;

	// Update annotations with new ids
	for (int energy = 0; energy < spectra.size(); energy++) {
		Spectrum::const_iterator it = spectra[energy].begin();
		for (int peak_idx = 0; it != spectra[energy].end(); ++it, peak_idx++) {
			for (unsigned int i = 0; i < it->annotations.size(); i++)
				spectra[energy].updateAnnotationId(peak_idx, i, id_lookup[it->annotations[i].first]);
		}
	}
}

void MolData::computeTransitionThetas(Param &param) {

	const unsigned int num_levels = param.getNumEnergyLevels();
	thetas.resize(num_levels);
	for (unsigned int energy = 0; energy < num_levels; energy++) {

		// Compute the theta value for each feature vector
		thetas[energy].resize(fg->getNumTransitions());
		for (unsigned int i = 0; i < fg->getNumTransitions(); i++)

			thetas[energy][i] = param.computeTheta(*(fg->getTransitionAtIdx(i)->getFeatureVector()), energy);
	}
}

void MolData::computeLogTransitionProbabilities() {

	log_probs.resize(thetas.size());
	for (unsigned int energy = 0; energy < thetas.size(); energy++) {

		// Will store an entry for all transitions, followed by all persistences
		log_probs[energy].resize(fg->getNumTransitions() + fg->getNumFragments());

		// Compute all the denominators
		std::vector<double> denom_cache(fg->getNumFragments());
		const tmap_t *from_id_map = fg->getFromIdTMap();
		for (unsigned int i = 0; i < fg->getNumFragments(); i++) {
			double denom = 1.0;
			for (auto trans_idx = (*from_id_map)[i].begin(); trans_idx != (*from_id_map)[i].end(); ++trans_idx) {
				denom += exp(thetas[energy][*trans_idx]);
			}
			denom_cache[i] = log(denom);
		}

		// Set the transition log probabilities
		for (unsigned int i = 0; i < fg->getNumTransitions(); i++) {
			const TransitionPtr t = fg->getTransitionAtIdx(i);
			// log(A/B) = log A - log B
			log_probs[energy][i]  = thetas[energy][i] - denom_cache[t->getFromId()];
		}

		// Set the persistence log probabilities
		int offset = fg->getNumTransitions();
		for (unsigned int i = 0; i < fg->getNumFragments(); i++) {
			if (fg->getFragmentAtIdx(i)->isIntermediate() && !cfg->allow_intermediate_peak)
				log_probs[energy][offset + i] = -100000000;
			else
				log_probs[energy][offset + i] = -denom_cache[i];
		}
	}
}

void MolData::readInSpectraFromFile(const std::string &peak_filename, bool readToPredicted) {

	// Set the spectra that we are writing to
	std::vector<Spectrum> *spec_dest = &spectra;
	if (readToPredicted) spec_dest = &predicted_spectra;

	// Clear any existing entries
	spec_dest->clear();
	spec_dest->resize(0);

	Spectrum *curr_spec;
	std::string line;
	std::ifstream ifs(peak_filename.c_str(), std::ifstream::in);
	if (!ifs) std::cerr << "Warning: Could not open file " << peak_filename << std::endl;

	// regex check
	std::regex peak_list_regex("^([0-9]+)([.][0-9]+)?([\\s\\t]{1,})([0-9]+)([.][0-9]+)?");

	int first = 1;
	while (ifs.good()) {
		getline(ifs, line);
		// 1.empty line separates peaks and annotations
		// 2.we need at least 4 char for peak: n space n newline
		if (line.length() < 2) break;
		// in case we are seen version string
		if (line[0] == '#') continue;
		// Check for the energy specifier - start a new spectrum if found
		// or start one anyway if there is no energy specifier
		if (line.substr(0, 3) == "low" || line.substr(0, 3) == "med" || line.substr(0, 4) == "high" ||
		    line.substr(0, 6) == "energy") {
			spec_dest->push_back(Spectrum());
			curr_spec = &(spec_dest->back());
			first     = 0;
			continue;
		} else if (first) {
			spec_dest->push_back(Spectrum());
			curr_spec = &(spec_dest->back());
		}
		first = 0;

		// Otherwise allocate peak to the current spectrum
		boost::trim(line);
		if (!std::regex_match(line, peak_list_regex)) {
			std::cerr << "Invalid peak " << line << std::endl;
			continue;
		}

		std::stringstream ss(line);
		double mass, intensity;
		ss >> mass >> intensity;
		curr_spec->push_back(Peak(mass, intensity));
	}
	ifs.close();

	for (auto &spec : *spec_dest) spec.normalizeAndSort();

	// once finished, copy spectrum to orig spectrum
	// because we already  sort and normalized
	orig_spectra.clear();
	for (const auto &spec : spectra) orig_spectra.push_back(spec);
}

void MolData::readInSpectraFromMSP(MspReader &msp, bool readToPredicted) {

	// Set the spectra that we are writing to
	std::vector<Spectrum> *spec_dest = &spectra;
	if (readToPredicted) spec_dest = &predicted_spectra;

	// Clear any existing entries
	spec_dest->clear();
	spec_dest->resize(0);

	// Copy the spectra in
	const std::vector<Spectrum> *spec = msp.fetchSpectrumForId(id.c_str());
	for (int energy = 0; energy < spec->size(); energy++) {
		spec_dest->push_back((*spec)[energy]);
		(*spec_dest)[energy].normalizeAndSort();
	}

	// once finished, copy specturm to orig
	// because we already  sort and normalized
	orig_spectra.clear();
	for (const auto &spectrum : spectra) orig_spectra.push_back(spectrum);
}

void MolData::cleanSpectra(double abs_tol, double ppm_tol) {

	for (auto &spectrum : spectra) spectrum.clean(abs_tol, ppm_tol);
}

std::string MolData::removePeaksWithNoFragment(double abs_tol, double ppm_tol) {

	std::vector<double> all_masses;
	getEnumerationSpectraMasses(all_masses);

	std::string msg;
	auto spectrum_idx = 1;
	for (auto &spectrum : spectra) {
		auto num_peaks   = spectrum.size();
		auto num_removed = spectrum.removePeaksWithNoFragment(all_masses, abs_tol, ppm_tol);
		msg += " Spectrum#" + std::to_string(spectrum_idx) + " removed: " + std::to_string(num_removed) + "/" +
		       std::to_string(num_peaks);
		spectrum_idx += 1;
	}
	return msg;
}

bool MolData::hasEmptySpectrum(int energy_level) const {

	bool result = false;
	if (energy_level < 0) {
		for (auto &spectrum : spectra) { result = (result || (spectrum.size() == 0)); }
	} else if (spectra.size() > energy_level) {
		result = spectra[energy_level].size() == 0;
	} else {
		result = true;
	}
	return result;
}

void MolData::computePredictedSpectra(Param &param, bool use_existing_thetas, int energy_level, int min_peaks,
                                      int max_peaks, double perc_thresh, double min_relative_intensity,
                                      int quantise_peaks_decimal_place, bool log_to_linear) {

	// Compute the transition probabilities using this parameter set
	if (!use_existing_thetas) computeTransitionThetas(param);
	computeLogTransitionProbabilities();

	// Generate and collect the peak results
	predicted_spectra.resize(cfg->spectrum_depths.size());
	if (energy_level == -1) {
		for (unsigned int energy = 0; energy < cfg->spectrum_depths.size(); energy++) {
			createSpeactraSingleEnergry(energy);
			if (log_to_linear) predicted_spectra[energy].convertToLinearScale();
		}
	} else {

		createSpeactraSingleEnergry(energy_level);
		if (log_to_linear) predicted_spectra[energy_level].convertToLinearScale();
	}

	postprocessPredictedSpectra(perc_thresh, min_peaks, max_peaks, min_relative_intensity,
	                            quantise_peaks_decimal_place);
}

void MolData::createSpeactraSingleEnergry(unsigned int energy_level) {
	predicted_spectra[energy_level].clear();

	config_t se_cfg;
	initSingleEnergyConfig(se_cfg, *cfg, energy_level);

	int depth = getFGHeight() > se_cfg.spectrum_depths[0] ? getFGHeight() : se_cfg.spectrum_depths[0];

	// Run forward inference
	std::vector<Message> msgs;
	Inference infer(this, &se_cfg);
	infer.runInferenceDownwardPass(msgs, depth, energy_level);

	// Extract the peaks from the relevant message
	int msg_depth = depth - 1; // se_cfg.spectrum_depths[0] - 1;
	Message *msg  = &(msgs[msg_depth]);
	if (cfg->include_isotopes)
		translatePeaksFromMsgToSpectraWithIsotopes(predicted_spectra[energy_level], msg);
	else
		translatePeaksFromMsgToSpectra(predicted_spectra[energy_level], msg);
}

void MolData::translatePeaksFromMsgToSpectra(Spectrum &out_spec, Message *msg) {

	// Create the peaks
	std::map<double, Peak> peak_probs;
	Message::const_iterator itt = msg->begin();
	for (; itt != msg->end(); ++itt) {
		double mass              = fg->getFragmentAtIdx(itt.index())->getMass();
		double intensity_contrib = exp(*itt);
		if (peak_probs.find(mass) != peak_probs.end())
			peak_probs[mass].intensity += intensity_contrib;
		else
			peak_probs[mass].intensity = intensity_contrib;
		peak_probs[mass].mass = mass;
		peak_probs[mass].annotations.push_back(annotation_t(itt.index(), intensity_contrib));
	}

	// Add the peaks to the spectra
	std::map<double, Peak>::iterator itm = peak_probs.begin();
	for (; itm != peak_probs.end(); ++itm) out_spec.push_back(itm->second);
}

void MolData::translatePeaksFromMsgToSpectraWithIsotopes(Spectrum &out_spec, Message *msg) {

	// Create the peaks
	std::map<double, Peak> peak_probs;
	Message::const_iterator itt = msg->begin();
	for (; itt != msg->end(); ++itt) {
		const Spectrum *isotope_spec = fg->getFragmentAtIdx(itt.index())->getIsotopeSpectrum();
		Spectrum::const_iterator itp = isotope_spec->begin();
		for (; itp != isotope_spec->end(); ++itp) {
			double intensity_contrib = 0.01 * itp->intensity * exp(*itt);
			double mass              = itp->mass;
			if (peak_probs.find(mass) != peak_probs.end())
				peak_probs[mass].intensity += intensity_contrib;
			else
				peak_probs[mass].intensity = intensity_contrib;
			peak_probs[mass].mass = mass;
			peak_probs[mass].annotations.push_back(annotation_t(itt.index(), intensity_contrib));
		}
	}

	// Add the peaks to the spectra
	std::map<double, Peak>::iterator itm = peak_probs.begin();
	for (; itm != peak_probs.end(); ++itm) out_spec.push_back(itm->second);
}

void MolData::writePredictedSpectraToFile(std::string &filename) {

	std::ofstream of;
	of.open(filename.c_str());
	if (!of.is_open()) {
		std::cout << "Warning: Trouble opening predicted spectrum file" << std::endl;
		return;
	}
	std::streambuf *buf = of.rdbuf();
	std::ostream out(buf);
	outputSpectra(out, "Predicted");
	of.close();
}

void MolData::writePredictedSpectraToMspFileStream(std::ostream &out) {
	for (unsigned int energy = 0; energy < getNumPredictedSpectra(); energy++) {
		const Spectrum *spec = getPredictedSpectrum(energy);
		spec->outputToMspStream(out, id, cfg->ionization_mode, energy, smiles_or_inchi, cfg->default_mz_decimal_place);
	}
}

void MolData::writePredictedSpectraToMgfFileStream(std::ostream &out) {
	double parent_ion_mw = getParentIonMass();

	for (unsigned int energy = 0; energy < getNumPredictedSpectra(); energy++) {
		const Spectrum *spec = getPredictedSpectrum(energy);
		spec->outputToMgfStream(out, id, cfg->ionization_mode, energy, parent_ion_mw, smiles_or_inchi,
		                        cfg->default_mz_decimal_place);
	}
}

void MolData::writeFullEnumerationSpectrumToFile(std::string &filename) {

	std::ofstream of;
	of.open(filename.c_str());
	if (!of.is_open()) {
		std::cout << "Warning: Trouble opening enumerated spectrum file" << std::endl;
		return;
	}
	std::streambuf *buf = of.rdbuf();
	std::ostream out(buf);
	std::vector<double> all_masses;
	getEnumerationSpectraMasses(all_masses);

	// All peaks of uniform intensity at each possible fragment mass
	for (unsigned int energy = 0; energy < spectra.size(); energy++) {

		out << "energy" << energy << std::endl;

		out << std::setprecision(10);
		double peak_height = 100.0 / (double)all_masses.size();
		for (int i = 0; i < all_masses.size(); i++) out << all_masses[i] << " " << peak_height << std::endl;
	}

	of.close();
}

void MolData::writeFullEnumerationSpectrumToMspFileStream(std::ostream &out) {

	std::vector<double> all_masses;
	getEnumerationSpectraMasses(all_masses);

	out << "ID: " << id << std::endl;
	out << "Num peaks: " << all_masses.size() << std::endl;
	out << std::setprecision(10);
	double peak_height = 100.0 / (double)all_masses.size();
	for (int i = 0; i < all_masses.size(); i++) out << all_masses[i] << " " << peak_height << std::endl;
	out << std::endl;
}

void MolData::getEnumerationSpectraMasses(std::vector<double> &output_masses) {

	// Get and sort the fragment masses
	unsigned int numf = fg->getNumFragments();
	std::vector<double> all_masses;
	if (fg->hasIsotopesIncluded()) {
		for (unsigned int i = 0; i < numf; i++) {

			const Fragment *f            = fg->getFragmentAtIdx(i);
			const Spectrum *isotope_spec = f->getIsotopeSpectrum();
			Spectrum::const_iterator itp = isotope_spec->begin();
			for (; itp != isotope_spec->end(); ++itp) all_masses.push_back(itp->mass);
		}
	} else {
		all_masses.resize(numf);
		for (unsigned int i = 0; i < numf; i++) {
			const Fragment *f = fg->getFragmentAtIdx(i);
			all_masses[i]     = f->getMass();
		}
	}
	std::sort(all_masses.begin(), all_masses.end());

	// Remove repeats
	output_masses.resize(all_masses.size());
	int i            = 0;
	double prev_mass = -1.0, tol = 1e-10;
	std::vector<double>::iterator it = all_masses.begin();
	for (; it != all_masses.end(); ++it) {
		if (fabs(*it - prev_mass) > tol) output_masses[i++] = *it;
		prev_mass = *it;
	}
	int num_unique = i;
	output_masses.resize(num_unique);
}

void MolData::outputSpectra(std::ostream &out, const char *spec_type, bool do_annotate, bool add_version) {

	std::vector<Spectrum> *spectra_to_output;
	if (std::string(spec_type) == "Predicted")
		spectra_to_output = &predicted_spectra;
	else if (std::string(spec_type) == "Experimental")
		spectra_to_output = &spectra;
	else
		std::cout << "Unknown spectrum type to output: " << spec_type << std::endl;

	if ((std::string(spec_type) == "Predicted") && add_version) {
		std::string spectra_str;
		switch (cfg->ionization_mode) {
		case (POSITIVE_ESI_IONIZATION_MODE): spectra_str = "ESI-MS/MS [M+H]+ Spectra"; break;

		case (NEGATIVE_ESI_IONIZATION_MODE): spectra_str = "ESI-MS/MS [M-H]- Spectra"; break;
		case (POSITIVE_EI_IONIZATION_MODE): spectra_str = "EI-MS Spectra"; break;
		default: break;
		}
		out << "#In-silico " << spectra_str << std::endl
		    << "#PREDICTED BY " << APP_STRING << " " << PROJECT_VER << std::endl;
		out << "#ID=" << this->getId() << std::endl;
		if (smiles_or_inchi.substr(0, 6) == "InChI=") {
			out << "#" << smiles_or_inchi << std::endl;
			out << "#InChiKey=" << RDKit::InchiToInchiKey(smiles_or_inchi) << std::endl;
		} else {
			RDKit::RWMol *rwmol;
			out << "#SMILES=" << smiles_or_inchi << std::endl;
			rwmol = RDKit::SmilesToMol(smiles_or_inchi);
			RDKit::ExtraInchiReturnValues rv;
			out << "#InChiKey=" << RDKit::InchiToInchiKey(RDKit::MolToInchi(*rwmol, rv)) << std::endl;
			out << "#Formula=" << RDKit::Descriptors::calcMolFormula(*rwmol) << std::endl;
			out << "#PMass=" << std::fixed << std::setprecision(5) << getParentIonMass() << std::endl;
			delete rwmol;
		}
	}

	auto it = spectra_to_output->begin();
	std::set<int> frag_ids;
	for (int energy = 0; it != spectra_to_output->end(); ++it, energy++) {
		out << "energy" << energy << std::endl;
		it->outputToStream(out, do_annotate, cfg->default_mz_decimal_place, true);
		it->getDisplayedFragmentIds(frag_ids);
	}

	if (do_annotate) {
		out << std::endl;
		writeFragmentsOnlyForIds(out, frag_ids);
	}
}

void MolData::postprocessPredictedSpectra(double perc_thresh, int min_peaks, int max_peaks,
                                          double min_relative_intensity_prec, int quantise_peaks_decimal_place) {

	for (auto &predicted_spectrum : predicted_spectra) {
		// predicted_spectrum.quantisePeaksByMass(10);
		if (quantise_peaks_decimal_place > -1) predicted_spectrum.quantisePeaksByMass(quantise_peaks_decimal_place);
		predicted_spectrum.postProcess(perc_thresh, min_peaks, max_peaks, min_relative_intensity_prec);
		predicted_spectrum.normalizeAndSort();
		predicted_spectrum.sortAndNormalizeAnnotations();
	}
}

void MolData::quantisePredictedSpectra(int num_dec_places) {

	for (auto &predicted_spectrum : predicted_spectra) { predicted_spectrum.quantisePeaksByMass(num_dec_places); }
}

void MolData::quantiseMeasuredSpectra(int num_dec_places) {

	for (auto &predicted_spectrum : predicted_spectra) { predicted_spectrum.quantisePeaksByMass(num_dec_places); }
}

void MolData::getSampledTransitionIdsWeightedRandomWalk(std::set<int> &selected_ids, int max_num_iter, int energy,
                                                        double explore_weight) {
	if (!hasEmptySpectrum(0) && hasComputedGraph())
		fg->getSampledTransitionIdsWeightedRandomWalk(selected_ids, max_num_iter, thetas[energy], explore_weight);
}

void MolData::getRandomSampledTransitions(std::set<int> &selected_ids, int max_selection) {
	if (!hasEmptySpectrum(0) && hasComputedGraph()) fg->getRandomSampledTransitions(selected_ids, max_selection);
}

void MolData::getSampledTransitionIdsRandomWalk(std::set<int> &selected_ids, int max_selection) {
	if (!hasEmptySpectrum(0) && hasComputedGraph()) fg->getSampledTransitionIdsRandomWalk(selected_ids, max_selection);
}

void MolData::getSampledTransitionIdUsingDiffMapBFS(std::set<int> &selected_ids,
                                                    std::set<unsigned int> &selected_weights) {
	if (!hasEmptySpectrum(0) && hasComputedGraph())
		fg->getSampledTransitionIdsDiffMapChildOnly(selected_ids, selected_weights);
}

void MolData::getSampledTransitionIdUsingDiffMapCA(std::set<int> &selected_ids,
                                                   std::set<unsigned int> &selected_weights) {
	if (!hasEmptySpectrum(0) && hasComputedGraph()) fg->getSampledTransitionIdsDiffMap(selected_ids, selected_weights);
}

// It is caller's response to compute predicted spectra
void MolData::getSelectedMasses(std::set<unsigned int> &selected_weights, int energry_level) {

	Comparator *cmp = new Dice(cfg->ppm_mass_tol, cfg->abs_mass_tol);
	std::vector<peak_pair_t> peak_pairs;
	cmp->getMatchingPeakPairsWithNoneMatchs(peak_pairs, &spectra[energry_level], &predicted_spectra[energry_level]);

	std::map<double, double, std::greater<double>> difference;
	for (const auto &peak_pair : peak_pairs) {
		double intensity_difference = std::fabs(peak_pair.first.intensity - peak_pair.second.intensity);
		// we are bias towards to peak should be there
		if (peak_pair.first.intensity > 0 && peak_pair.second.intensity == 0) intensity_difference += 100.0;
		if (intensity_difference > cfg->ga_diff_sampling_difference) {
			double peak_mass = peak_pair.second.mass;
			difference.insert(std::pair<double, double>(intensity_difference, peak_mass));
		}
	}
	delete (cmp);

	for (const auto &diff : difference) {
		if (selected_weights.size() >= cfg->ga_diff_sampling_peak_num) break;
		// we only need 5 digitals after decimal
		selected_weights.insert((unsigned int)std::round(diff.second * WEIGHT_SELECTION_SCALER));
	}
}

void MolData::computeMergedPrediction() {
	if (m_merged_predicted_spectra == nullptr) m_merged_predicted_spectra = new Spectrum();

	std::map<double, double> peaks_map;
	for (auto &spectrum : predicted_spectra) {
		for (auto &peak : *spectrum.getPeaks()) {
			if (peaks_map.find(peak.mass) == peaks_map.end()) peaks_map[peak.mass] = 0.0;
			peaks_map[peak.mass] += peak.intensity;
		}
	}
	for (const auto &peak_data : peaks_map) {
		Peak peak(peak_data.first, peak_data.second);
		m_merged_predicted_spectra->push_back(peak);
	}
	m_merged_predicted_spectra->normalizeAndSort();
}

double MolData::getParentIonMass() const {
	romol_ptr_t mol    = createMolPtr(smiles_or_inchi.c_str());
	double parent_mass = getMonoIsotopicMass(mol);
	switch (cfg->ionization_mode) {
	case (POSITIVE_ESI_IONIZATION_MODE): parent_mass += 1.007276; break;
	case (NEGATIVE_ESI_IONIZATION_MODE): parent_mass -= 1.007276; break;
	case (POSITIVE_EI_IONIZATION_MODE): parent_mass -= MASS_ELECTRON; break;
	default: break;
	}
	return parent_mass;
}

MolData::~MolData() {

	delete fg;
	delete ev_fg;
	delete m_merged_predicted_spectra;
}
