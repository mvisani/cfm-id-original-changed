/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM_NN.cpp
#
# Description: 	Class to apply Expectation Maximization algorithm to derive
#				model parameters when using a neural net for theta.
#				- all identical to linear model, except params are NNParam and gradient
#				  computation is different.
#					E-step: IPFP or equivalent.
#					M-step: Gradient Ascent
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "EmNNModel.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cstddef>
#include <iostream>

EmNNModel::EmNNModel(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename,
                     std::string initial_params_filename)
    : EmModel(a_cfg, an_fc, a_status_filename, initial_params_filename) {
	int num_energies_to_include = cfg->map_d_to_energy.back() + 1;
	if (!initial_params_provided) {
		nn_param = boost::shared_ptr<NNParam>(new NNParam(
		    fc->getFeatureNames(), num_energies_to_include, a_cfg->theta_nn_hlayer_num_nodes,
		    a_cfg->theta_nn_layer_act_func_ids, a_cfg->nn_layer_dropout_probs, a_cfg->nn_layer_is_frozen_flags));
		std::cout << "EM_NN: No initial params provided" << std::endl;
	} else {
		nn_param = boost::shared_ptr<NNParam>(new NNParam(initial_params_filename));
		while (nn_param->getNumEnergyLevels() < num_energies_to_include) nn_param->appendRepeatedPrevEnergyParams();
		std::string msg = "EM_NN: Initial params provided from " + initial_params_filename;
		std::cout << msg << std::endl;
	}
	this->param = nn_param;
}

// These functions are the same as for EM, but use NNParam, so need to be redefined here
void EmNNModel::computeThetas(MolData *moldata) { moldata->computeTransitionThetas(*nn_param); }

void EmNNModel::writeParamsToFile(std::string &filename) { nn_param->saveToFile(filename); }

// Gradient Computation using Backpropagation
int EmNNModel::computeAndAccumulateGradient(double *grads, int mol_idx, MolData &mol_data, suft_counts_t &suft,
                                            int sampling_method, unsigned int energy) {
	unsigned int num_transitions    = mol_data.getNumTransitions();
	unsigned int num_fragments      = mol_data.getNumFragments();
	unsigned int offset             = num_transitions;
	unsigned int layer2_offset      = nn_param->getSecondLayerWeightOffset();
	unsigned int weights_per_energy = nn_param->getNumWeightsPerEnergyLevel();

	size_t num_used_transitions = 0;
	if (!mol_data.hasComputedGraph()) return static_cast<int>(num_used_transitions);

	suft_t *suft_values = &(suft.values[mol_idx]);

	unsigned int grad_offset = energy * nn_param->getNumWeightsPerEnergyLevel();
	unsigned int suft_offset = energy * (num_transitions + num_fragments);

	std::set<int> selected_trans_id;
	if (sampling_method != USE_NO_SAMPLING) {
		getSubSampledTransitions(mol_data, sampling_method, energy, selected_trans_id);
		num_used_transitions = selected_trans_id.size();
	}

	// Iterate over from_id (i)
	const tmap_t *from_map = mol_data.getFromIdTMap();
	for (int from_idx = 0; from_idx < num_fragments; from_idx++) {
		const std::vector<int> *from_id_map = &(*from_map)[from_idx];
		std::vector<int> sampled_trans_id;

		if (sampling_method != USE_NO_SAMPLING) {
			for (const auto &trans_id : (*from_map)[from_idx])
				if (selected_trans_id.find(trans_id) != selected_trans_id.end()) sampled_trans_id.push_back(trans_id);
			from_id_map = &sampled_trans_id;
		}

		if (from_id_map->empty()) continue;

		unsigned int num_trans_from_id = from_id_map->size();
		std::vector<azd_vals_t> a_values(num_trans_from_id);
		std::vector<azd_vals_t> z_values(num_trans_from_id);

		// Compute the forward values, storing intermediate a and z values, and the combined denom of the rho term, and
		// Q
		double denom = 1.0;
		auto it      = from_id_map->begin();
		std::vector<const FeatureVector *> fvs(num_trans_from_id);
		std::vector<double> nu_terms(num_trans_from_id + 1);
		for (int idx = 0; it != from_id_map->end(); ++it, idx++) {
			fvs[idx]     = mol_data.getFeatureVectorForIdx(*it);
			// use nn_param->compute theta in forward pass mode
			// which uses Inverted Dropout
			// do not use drop outs during used idxs collection
			// Otherwise this will cause segfault
			// You have been warned
			double theta = nn_param->computeTheta(*fvs[idx], energy, z_values[idx], a_values[idx], false, true);
			denom += exp(theta);
			nu_terms[idx] = (*suft_values)[*it + suft_offset];
		}
		nu_terms[num_trans_from_id] = (*suft_values)[offset + from_idx + suft_offset];

		// Compute the deltas
		std::vector<azd_vals_t> deltasA, deltasB;
		nn_param->computeDeltas(deltasA, deltasB, z_values, a_values, denom, energy);
		// Compute the unweighted gradients
		std::vector<std::vector<float>> unweighted_grads;
		std::set<unsigned int> from_id_used_idxs;
		nn_param->computeUnweightedGradients(unweighted_grads, from_id_used_idxs, fvs, deltasA, deltasB, a_values);

		// Accumulate the weighted gradients
		std::set<unsigned int>::iterator sit;
		it = from_id_map->begin();
		for (int idx = 0; idx <= num_trans_from_id; idx++) {
			// First layer
			double nu = nu_terms[idx];
			for (sit = from_id_used_idxs.begin(); sit != from_id_used_idxs.end(); ++sit)
				*(grads + grad_offset + *sit) += nu * unweighted_grads[idx][*sit];
			// Remaining layers
			for (int i = layer2_offset; i < weights_per_energy; i++)
				*(grads + grad_offset + i) += nu * unweighted_grads[idx][i];
		}
	}
	return num_used_transitions;
}

void EmNNModel::collectUsedIdx(MolData &mol_data, std::set<unsigned int> &used_idxs, unsigned int energy) {

	if (!mol_data.hasComputedGraph()) return;

	unsigned int num_fragments = mol_data.getNumFragments();
	unsigned int grad_offset   = energy * nn_param->getNumWeightsPerEnergyLevel();

	// Iterate over from_id (i)
	for (auto frag_trans_map = mol_data.getFromIdTMap()->begin(); frag_trans_map != mol_data.getFromIdTMap()->end();
	     ++frag_trans_map) {
		for (auto trans_id : *frag_trans_map) {
			const FeatureVector *fv  = mol_data.getFeatureVectorForIdx(trans_id);
			unsigned int feature_len = fv->getTotalLength();
			for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it)
				nn_param->collectUsedIdx(used_idxs, feature_len, grad_offset, *fv_it, energy);
		}
	}

	int weights_per_energy     = nn_param->getNumWeightsPerEnergyLevel();
	unsigned int layer2_offset = nn_param->getSecondLayerWeightOffset();

	// Add used indexes for the subsequent layers (all of them)
	for (int i = layer2_offset; i < weights_per_energy; i++) used_idxs.insert(grad_offset + i);
}

double EmNNModel::computeLogLikelihoodLoss(int molidx, MolData &moldata, suft_counts_t &suft, unsigned int energy) {

	double Q                     = 0.0;
	unsigned int num_transitions = moldata.getNumTransitions();
	unsigned int num_fragments   = moldata.getNumFragments();

	int offset = num_transitions;

	if (!moldata.hasComputedGraph()) return Q;

	suft_t *suft_values = &(suft.values[molidx]);

	unsigned int suft_offset = energy * (num_transitions + num_fragments);

	// Iterate over from_id (i)
	const tmap_t *from_map = moldata.getFromIdTMap();
	for (int from_idx = 0; from_idx < num_fragments; from_idx++) {

		const std::vector<int> *from_id_map = &(*from_map)[from_idx];
		unsigned int num_trans_from_id      = from_id_map->size();

		// Compute the forward values, and the combined denom of the rho term, and Q
		double denom = 1.0, nu_sum = 0.0;
		auto it = from_id_map->begin();
		std::vector<const FeatureVector *> fvs(num_trans_from_id);
		std::vector<double> nu_terms(num_trans_from_id + 1);
		for (int idx = 0; it != from_id_map->end(); ++it, idx++) {
			fvs[idx]     = moldata.getFeatureVectorForIdx(*it);
			double theta = nn_param->computeTheta(*fvs[idx], energy);
			denom += exp(theta);
			nu_terms[idx] = (*suft_values)[*it + suft_offset];
			nu_sum += nu_terms[idx];
			Q += nu_terms[idx] * theta;
		}
		nu_terms[num_trans_from_id] = (*suft_values)[offset + from_idx + suft_offset];
		nu_sum += nu_terms[num_trans_from_id];
		Q -= log(denom) * nu_sum;
	}
	return Q;
}

void EmNNModel::updateGradientForRegularizationTerm(double *grads, unsigned int energy) {

	auto it = this->used_idxs.begin();
	for (; it != this->used_idxs.end(); ++it) {
		double weight = nn_param->getWeightAtIdx(*it);
		*(grads + *it) -= cfg->lambda * weight;
	}

	// Remove part of the Bias terms (regularize the bias terms 100 times less than the other weights!)
	unsigned int weights_per_energy = nn_param->getNumWeightsPerEnergyLevel();
	std::vector<unsigned int> bias_indexes;
	nn_param->getBiasIndexes(bias_indexes);
	int offset = energy * weights_per_energy;
	for (auto it = bias_indexes.begin(); it != bias_indexes.end(); ++it) {
		double bias = nn_param->getWeightAtIdx(offset + *it);
		*(grads + offset + *it) += 0.99 * cfg->lambda * bias;
	}
}

double EmNNModel::getRegularizationTerm(unsigned int energy) {

	double req_term = 0.0;
	auto it         = this->used_idxs.begin();
	for (; it != this->used_idxs.end(); ++it) {
		double weight = nn_param->getWeightAtIdx(*it);
		req_term -= 0.5 * cfg->lambda * weight * weight;
	}

	// Remove part of the Bias terms (regularize the bias terms 100 times less than the other weights!)
	unsigned int weights_per_energy = nn_param->getNumWeightsPerEnergyLevel();
	std::vector<unsigned int> bias_indexes;
	nn_param->getBiasIndexes(bias_indexes);
	int offset = energy * weights_per_energy;
	for (auto it = bias_indexes.begin(); it != bias_indexes.end(); ++it) {
		double bias = nn_param->getWeightAtIdx(offset + *it);
		req_term += 0.99 * 0.5 * cfg->lambda * bias * bias;
	}
	return req_term;
}
