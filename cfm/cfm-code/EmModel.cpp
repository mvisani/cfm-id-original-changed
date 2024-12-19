/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.cpp
#
# Description: 	Class to apply Expectation Maximization algorithm to derive
#				model parameters.
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
#include "Config.h"
#include "omp.h"
#include <boost/filesystem.hpp>
#include <cfloat>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <ostream>
#include <tuple>

#include "Comparators.h"
#include "EmModel.h"
#include "Inference.h"
#include "Param.h"

// constructor
EmModel::EmModel(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename,
                 const std::string &initial_params_filename) {
	this->cfg             = a_cfg;
	this->fc              = an_fc;
	this->status_filename = a_status_filename;

	int num_energies_to_include = cfg->map_d_to_energy.back() + 1;
	if (initial_params_filename.empty()) {
		param                   = boost::shared_ptr<Param>(new Param(fc->getFeatureNames(), num_energies_to_include));
		initial_params_provided = false;
		std::cout << "EM: No initial params provided" << std::endl;
	} else {
		param = boost::shared_ptr<Param>(new Param(initial_params_filename));
		while (param->getNumEnergyLevels() < num_energies_to_include) { param->appendRepeatedPrevEnergyParams(); }
		initial_params_provided = true;
		std::string msg         = "EM: Initial params provided from " + initial_params_filename;
		std::cout << msg << std::endl;
	}
	sparse_params = true;
	start_time    = std::chrono::system_clock::now();
}

void EmModel::computeThetas(MolData *moldata) { moldata->computeTransitionThetas(*param); }

std::tuple<double, double, double, double> EmModel::computeMetrics(const Spectrum *p, const Spectrum *q) {
	double dice_out, dp_out, precision_out, recall_out;

	Comparator *dice_cmp       = new Dice(cfg->ppm_mass_tol, cfg->abs_mass_tol);
	Comparator *dotproduct_cmp = new DotProduct(cfg->ppm_mass_tol, cfg->abs_mass_tol);
	Comparator *p_cmp          = new Precision(cfg->ppm_mass_tol, cfg->abs_mass_tol);
	Comparator *r_cmp          = new Recall(cfg->ppm_mass_tol, cfg->abs_mass_tol);

	dice_out      = dice_cmp->computeScore(p, q);
	dp_out        = dotproduct_cmp->computeScore(p, q);
	precision_out = p_cmp->computeScore(p, q);
	recall_out    = r_cmp->computeScore(p, q);

	delete dice_cmp;
	delete dotproduct_cmp;
	delete p_cmp;
	delete r_cmp;
	return std::make_tuple(dice_out, dp_out, precision_out, recall_out);
}

void EmModel::initSuft(suft_counts_t &suft, std::vector<MolData> &data) {

	// Resize the suft structure for each molecule
	auto num_mols = data.size();
	suft.values.resize(num_mols);
	for (unsigned int i = 0; i < num_mols; i++) {
		// const FragmentGraph *fg = data[i].getFragmentGraph();
		size_t len         = data[i].getNumTransitions() + data[i].getNumFragments();
		size_t num_spectra = data[i].getNumSpectra();
		suft.values[i].resize(len * num_spectra);
	}
}

void EmModel::recordSufficientStatistics(suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs,
                                         unsigned int energy) {

	int depth = cfg->model_depth > moldata->getFGHeight() ? cfg->model_depth : moldata->getFGHeight();
	unsigned int num_transitions = moldata->getNumTransitions();
	unsigned int num_fragments   = moldata->getNumFragments();
	unsigned int len_offset      = num_transitions + num_fragments;

	// Accumulate the Sufficient Statistics
	for (unsigned int i = 0; i < num_transitions; i++) {

		const TransitionPtr t = moldata->getTransitionAtIdx(i);

		double belief = 0.0;
		// int energy = cfg->map_d_to_energy[0];
		if (t->getFromId() == 0) // main ion is always id = 0
			belief += exp(beliefs->tn[i][0]);

		for (unsigned int d = 1; d < depth; d++) { belief += exp(beliefs->tn[i][d]); }
		suft.values[molidx][i + energy * len_offset] = belief;
	}

	// Accumulate the persistence terms
	unsigned int offset = num_transitions;
	for (unsigned int i = 0; i < num_fragments; i++) {

		double belief = 0.0;

		if (i == 0) // main ion is always id = 0
			belief += exp(beliefs->ps[i][0]);
		for (unsigned int d = 1; d < depth; d++) { belief += exp(beliefs->ps[i][d]); }
		suft.values[molidx][i + offset + energy * len_offset] = belief;
	}
}

void EmModel::collectUsedIdx(MolData &mol_data, std::set<unsigned int> &used_idxs, unsigned int energy) {

	if (!mol_data.hasComputedGraph()) return;

	unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();

	// Iterate over from_id (i)
	for (auto frag_trans_map = mol_data.getFromIdTMap()->begin(); frag_trans_map != mol_data.getFromIdTMap()->end();
	     ++frag_trans_map) {
		for (auto trans_id : *frag_trans_map) {
			const FeatureVector *fv = mol_data.getFeatureVectorForIdx(trans_id);
			for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it)
				used_idxs.insert(*fv_it + grad_offset);
		}
	}
}

void EmModel::zeroUnusedParams() {
	unsigned int i;
	for (i = 0; i < param->getNumWeights(); i++) {
		if (used_idxs.find(i) == used_idxs.end()) param->setWeightAtIdx(0.0, i);
	}
}

double EmModel::updateParametersGradientAscent(std::vector<MolData> &data, suft_counts_t &suft, double learning_rate,
                                               int sampling_method, unsigned int energy) {
	// DBL_MIN is the smallest positive double
	// -DBL_MAX is the smallest negative double
	double loss = 0.0, prev_loss = -DBL_MAX, prev_best_loss = -DBL_MAX;

	std::vector<double> grads(this->param->getNumWeights(), 0.0);
	Solver *solver = nullptr;
	solver         = getSolver(cfg->ga_method, learning_rate);

	if (this->used_idxs.empty()) {
		auto before = std::chrono::system_clock::now();
		if (cfg->collected_all_used_idx) {
			auto grad_offset = energy * param->getNumWeightsPerEnergyLevel();
			for (int i = 0; i < param->getNumWeightsPerEnergyLevel(); ++i) { this->used_idxs.insert(grad_offset + i); }
		}

#pragma omp parallel for
		for (size_t molidx = 0; molidx < data.size(); ++molidx) {
			MolData &mol = data[molidx];
			if (mol.getGroup() != validation_group) {
				if (!cfg->collected_all_used_idx) {
					std::set<unsigned int> local_used_idxs;
					collectUsedIdx(mol, local_used_idxs, energy);

// Merge the local set of used indices into the shared set
#pragma omp critical
					{ used_idxs.insert(local_used_idxs.begin(), local_used_idxs.end()); }
				}
				computeThetas(&mol);
			}
		}
		zeroUnusedParams();
		auto after = std::chrono::system_clock::now();
		std::cout << "[M-Step][T+" << getTimeDifferenceStr(start_time, after)
		          << "s]Collect Used Index, Time Used: " << getTimeDifferenceStr(before, after) + "s" << std::endl;
	}

	int iter                 = 0;
	int max_iteration        = cfg->ga_max_iterations;
	int ga_no_progress_count = 0;
	std::vector<float> current_best_weight;
	while (iter++ < max_iteration && ga_no_progress_count <= cfg->ga_no_progress_count) {
		auto before = std::chrono::system_clock::now();
		if (iter > 1) { prev_loss = loss; }
		if (USE_NO_DECAY != cfg->ga_decay_method) {
			learning_rate = gaGetUpdatedLearningRate(learning_rate, iter);
			solver->setLearningRate(learning_rate);
		}
		// Select molecules to include in gradient mini-batch.
		std::vector<int> minibatch_flags(data.size());
		int num_batch = cfg->ga_minibatch_nth_size;
		setMiniBatchFlags(minibatch_flags, num_batch);

		// Compute the gradient
		std::clock_t c_start = clock();

		for (auto batch_idx = 0; batch_idx < num_batch; ++batch_idx) {
			double num_trans = 0;
			int molidx       = 0;

			std::fill(grads.begin(), grads.end(), 0.0);

#pragma omp parallel for reduction(+ : num_trans)
			for (int molidx = 0; molidx < data.size(); ++molidx) {
				auto &mol_it = data[molidx];

				if (minibatch_flags[molidx] == batch_idx && mol_it.getGroup() != validation_group) {
					// so now it should not crash anymore
					if (!mol_it.hasComputedGraph()) continue;

					// Store the local gradient for the current thread
					std::vector<double> local_grads(grads.size(), 0.0);
					num_trans +=
					    computeAndAccumulateGradient(&local_grads[0], molidx, mol_it, suft, sampling_method, energy);

#pragma omp critical
					{
						for (size_t i = 0; i < grads.size(); ++i) { grads[i] += local_grads[i]; }
					}
				}
			}
			// update L2 only if lambda > 0
			if (cfg->lambda > 0.0) updateGradientForRegularizationTerm(&grads[0], energy);
			if (num_trans > 0)
				for (auto &grad : grads) grad /= num_trans;
			solver->adjustWeights(grads, this->used_idxs, param);
		}

		// End of epoch

		std::clock_t c_end = clock();
		std::string cpu_usage_string;

		// compute loss
		loss = computeAndSyncLoss(data, suft, energy);

		auto after            = std::chrono::system_clock::now();
		auto loss_change_rate = 1.0 - loss / prev_best_loss;
		std::cout << iter << ".[T+" << getTimeDifferenceStr(start_time, after) << "s]" << "Loss=" << loss
		          << " Prev=" << prev_loss << " Best=" << prev_best_loss << " Change Rate=" << loss_change_rate
		          << " Time Escaped=" << getTimeDifferenceStr(before, after) << "s Learning Rate=" << learning_rate;

		// let us roll Dropouts
		param->rollDropouts();

		// compute if we need stop
		if (loss_change_rate < cfg->ga_converge_thresh)
			ga_no_progress_count++;
		else
			ga_no_progress_count = 0;

		// we are doing ga, save best weight some where and check if we have getting better
		if (prev_best_loss < loss) {
			std::cout << " [Best Param]";
			current_best_weight = *param->getWeightsPtr();
			prev_best_loss      = loss;
		}
		std::cout << std::endl;
	}

	if (iter == cfg->ga_max_iterations)
		std::cout << "[M-Step]Gradient ascent did not converge" << std::endl;
	else
		std::cout << "[M-Step]Gradient ascent converged after " << iter - 1 << " iterations" << std::endl;
	delete solver;

	param->setWeights(current_best_weight);
	loss = prev_best_loss;

	return loss;
}

double EmModel::trainModel(std::vector<MolData> &molDataSet, int group, std::string &out_param_filename,
                           int energy_level) {
	unused_zeroed = 0;
	int iter      = 0;

	validation_group                    = group;
	// Write the initialised params to file (we may get want to reload and use
	// with saved suft state, even before updating)
	std::string init_out_param_filename = out_param_filename + "_init";
	writeParamsToFile(init_out_param_filename);

	// EM
	iter = 0;
	double loss;
	double prev_loss = -DBL_MAX;
	double best_loss = -DBL_MAX;

	// make of copy of learing rate
	// so we can share the save lr var over all em iterations
	// init some flags
	float learning_rate      = cfg->starting_step_size;
	int sampling_method      = cfg->ga_sampling_method;
	int em_no_progress_count = 0;

	int num_training_mols = 0, num_val_mols = 0;

	double max_pob_dice = 0.0, max_pob_dp = 0.0, max_pob_precision = 0.0, max_pob_recall = 0.0;

#pragma omp parallel for reduction(+ : max_pob_dice, max_pob_dp, max_pob_precision, max_pob_recall, num_training_mols, \
                                       num_val_mols)
	for (auto &mol_data : molDataSet) {
#pragma omp critial
		if (mol_data.hasEmptySpectrum(energy_level) && mol_data.getGroup() != validation_group)
			std::cout << "Warning: No peaks with explanatory fragment found for " << mol_data.getId()
			          << ", ignoring this input molecule." << std::endl;

		// compute reachable mectirs with curret fgs
		if (mol_data.getGroup() != validation_group) {
			auto [internal_prob_dice, internal_pob_db, internal_precision, internal_recall] =
			    computeMetrics(mol_data.getOrigSpectrum(energy_level), mol_data.getSpectrum(energy_level));
			max_pob_dice += internal_prob_dice;
			max_pob_dp += internal_pob_db;
			max_pob_precision += internal_precision;
			max_pob_recall += internal_recall;
			num_training_mols += 1;

		} else {
			num_val_mols += 1;
		}
		if (cfg->use_log_scale_peak) mol_data.convertSpectraToLogScale();
	}

	std::string pre_train_str;
	pre_train_str += "[Pre-Train][Max Possible Metrics for Current Fragmentation Setting]\n";
	pre_train_str += "Dice=" + std::to_string(max_pob_dice / num_training_mols) +
	                 " DotProduct=" + std::to_string(max_pob_dp / num_training_mols) +
	                 " Precision=" + std::to_string(max_pob_precision / num_training_mols / 100.0f) +
	                 " Recall=" + std::to_string(max_pob_recall / num_training_mols / 100.0f);

	std::cout << pre_train_str << std::endl;

	while (iter < cfg->em_max_iterations) {
		std::string iter_out_param_filename = out_param_filename + "_" + std::to_string(iter);
		std::string msg                     = "[Training]EM Iteration " + std::to_string(iter);
		std::cout << msg << std::endl;

		// Reset sufficient counts
		suft_counts_t suft;
		initSuft(suft, molDataSet);

		auto before = std::chrono::system_clock::now();

#pragma omp parallel for
		for (int molidx = 0; molidx < molDataSet.size(); ++molidx) {
			auto &mol = molDataSet[molidx];
			if (!mol.hasComputedGraph()) { continue; }
			if (mol.hasEmptySpectrum(energy_level)) { continue; }
			if (mol.getGroup() == validation_group && cfg->disable_cross_val_metrics) { continue; }
			computeThetas(&mol);
			mol.computeLogTransitionProbabilities();

			// Apply the peak evidence, compute the beliefs and record the sufficient
			// statistics
			beliefs_t beliefs;
			Inference infer(&mol, cfg);
			infer.calculateBeliefs(beliefs, energy_level);

#pragma omp critical // since we are updating shared suft we need to make sure we are not writing to the same memory
                     // location
			recordSufficientStatistics(suft, molidx, &mol, &beliefs, energy_level);
		}

		auto after = std::chrono::system_clock::now();
		std::string estep_time_msg =
		    "[E-Step][T+" + getTimeDifferenceStr(start_time, after) +
		    "s]Completed E-step processing: Time Elapsed = " + getTimeDifferenceStr(before, after) + " s";
		std::cout << estep_time_msg << std::endl;
		before = std::chrono::system_clock::now();
		std::cout << "[M-Step]Staring Learning Rate=" << learning_rate << std::endl;
		loss  = updateParametersGradientAscent(molDataSet, suft, learning_rate, sampling_method, energy_level);
		after = std::chrono::system_clock::now();
		std::string param_update_time_msg =
		    "[M-Step][T+" + std::to_string(getTimeDifference(start_time, after)) +
		    "s]Completed M-step update: Time Elapsed = " + getTimeDifferenceStr(before, after);
		std::cout << param_update_time_msg << std::endl;

		before = std::chrono::system_clock::now();

		// validation Q value
		double val_q = 0.0;

		int molidx        = 0;
		double train_dice = 0.0, train_dp = 0.0, train_precision = 0.0, train_recall = 0.0;
		double val_dice = 0.0, val_dp = 0.0, val_precision = 0.0, val_recall = 0.0;

		double pruned_train_dice = 0.0, pruned_train_dp = 0.0, pruned_train_precision = 0.0, pruned_train_recall = 0.0;

#pragma omp parallel for reduction(+ : val_q, train_dice, train_dp, train_precision, train_recall, val_dice, val_dp,   \
                                       val_precision, val_recall, pruned_train_dice, pruned_train_dp,                  \
                                       pruned_train_precision, pruned_train_recall)
		for (int molidx = 0; molidx < molDataSet.size(); ++molidx) {
			auto &mol_it = molDataSet[molidx];
			// run preduiction
			mol_it.computePredictedSpectra(*param, false, energy_level, cfg->default_predicted_peak_min,
			                               cfg->default_predicted_peak_max, cfg->default_postprocessing_energy,
			                               cfg->default_predicted_min_intensity, cfg->default_mz_decimal_place,
			                               cfg->use_log_scale_peak);

			if (mol_it.getGroup() == validation_group && !cfg->disable_cross_val_metrics) {
				// num_val_mols++;
				val_q += computeLogLikelihoodLoss(molidx, mol_it, suft, energy_level);

				auto [val_dice_internal, val_dp_internal, val_precision_internal, val_recall_internal] =
				    computeMetrics(mol_it.getOrigSpectrum(energy_level), mol_it.getPredictedSpectrum(energy_level));
				val_dice += val_dice_internal;
				val_dp += val_dp_internal;
				val_precision += val_precision_internal;
				val_recall += val_recall_internal;

			} else if (mol_it.getGroup() != validation_group && !cfg->disable_training_metrics) {
				auto [train_dice_internal, train_dp_internal, train_precision_internal, train_recall_internal] =
				    computeMetrics(mol_it.getOrigSpectrum(energy_level), mol_it.getPredictedSpectrum(energy_level));

				train_dice += train_dice_internal;
				train_dp += train_dp_internal;
				train_precision += train_precision_internal;
				train_recall += train_recall_internal;
				auto [pruned_train_dice_internal, pruned_train_dp_internal, pruned_train_precision_internal,
				      pruned_train_recall_internal] =
				    computeMetrics(mol_it.getSpectrum(energy_level), mol_it.getPredictedSpectrum(energy_level));
				pruned_train_dice += pruned_train_dice_internal;
				pruned_train_dp += pruned_train_dp_internal;
				pruned_train_precision += pruned_train_precision_internal;
				pruned_train_recall += pruned_train_recall_internal;
			}
		}

		after                  = std::chrono::system_clock::now();
		std::string q_time_msg = "[M-Step][T+" + getTimeDifferenceStr(start_time, after) +
		                         "s]Finished loss compute: Time Elapsed = " + getTimeDifferenceStr(before, after) +
		                         " seconds";
		std::cout << q_time_msg << std::endl;

		// Check for convergence
		// double loss_ratio = fabs((loss - prev_loss) / loss);
		auto loss_change_rate = 1.0 - loss / prev_loss;

		std::string qdif_str;
		qdif_str += "[M-Step][Traning Loss]    ";
		qdif_str += "Total=" + std::to_string(loss) + " Mean=" + std::to_string(loss / num_training_mols);

		if (prev_loss != -DBL_MAX)
			qdif_str += " Change Ratio= " + std::to_string(loss_change_rate) + " Prev=" + std::to_string(prev_loss);
		if (best_loss != -DBL_MAX) qdif_str += " PrevBest=" + std::to_string(best_loss);

		if (!cfg->disable_cross_val_metrics) {
			qdif_str += "\n[M-Step][Validation Loss (Without L2 Reg)] ";
			qdif_str += "Total=" + std::to_string(val_q) + " Mean=" + std::to_string(val_q / num_val_mols);
		}

		if (!cfg->disable_training_metrics) {
			qdif_str += "\n[M-Step][Training Metric]          ";
			qdif_str += "Dice=" + std::to_string(train_dice / num_training_mols) +
			            " DotProduct=" + std::to_string(train_dp / num_training_mols) +
			            " Precision=" + std::to_string(train_precision / num_training_mols / 100.0f) +
			            " Recall=" + std::to_string(train_recall / num_training_mols / 100.0f);

			qdif_str += "\n[M-Step][Training Metric (Pruned)] ";
			qdif_str += "Dice=" + std::to_string(pruned_train_dice / num_training_mols) +
			            " DotProduct=" + std::to_string(pruned_train_dp / num_training_mols) +
			            " Precision=" + std::to_string(pruned_train_precision / num_training_mols / 100.0f) +
			            " Recall=" + std::to_string(pruned_train_recall / num_training_mols / 100.0f);
		}

		if (!cfg->disable_cross_val_metrics) {
			qdif_str += "\n[M-Step][Validation Metric]        ";
			qdif_str += "Dice=" + std::to_string(val_dice / num_val_mols) +
			            " DotProduct=" + std::to_string(val_dp / num_val_mols) +
			            " Precision=" + std::to_string(val_precision / num_val_mols / 100.0f) +
			            " Recall=" += std::to_string(val_recall / num_val_mols / 100.0f);
		}
		std::cout << qdif_str << std::endl;

		if (best_loss < loss) {
			best_loss = loss;
			std::string progress_str =
			    "[M-Step] Found Better Q: " + std::to_string(best_loss) + " Write to File: " + out_param_filename;
			std::cout << progress_str << std::endl;
			writeParamsToFile(out_param_filename);
		}
		writeParamsToFile(iter_out_param_filename);
		updateEmTrainingParams(loss_change_rate, learning_rate, sampling_method, em_no_progress_count);
		// std::cerr << em_no_progress_count << std::endl;
		if (em_no_progress_count >= cfg->em_no_progress_count) {

			std::cout << "EM Stopped after " + std::to_string(em_no_progress_count) + " No Progress Iterations"
			          << std::endl;
			;
			std::cout << "EM Converged after " + std::to_string(iter) + " iterations" << std::endl;

			// time to stop
			// before stop, let us load best model
			param->readFromFile(out_param_filename);
			break;
		}
		prev_loss = loss;
		iter++;
	}
	if (iter >= cfg->em_max_iterations) {
		std::cout << "Warning: EM still converging after " + std::to_string(iter) + " iterations." << std::endl;
	}
	return best_loss;
}

float EmModel::getTimeDifference(const std::chrono::system_clock::time_point &before,
                                 const std::chrono::system_clock::time_point &after) {
	auto count = std::chrono::duration_cast<std::chrono::milliseconds>(after - before).count();
	return count / 1000.0f;
}

std::string EmModel::getTimeDifferenceStr(const std::chrono::system_clock::time_point &before,
                                          const std::chrono::system_clock::time_point &after) {
	float value = getTimeDifference(before, after);
	std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << value;
	return stream.str();
}

void EmModel::updateEmTrainingParams(double loss_change_rate, float &learning_rate, int &sampling_method,
                                     int &count_no_progress) const {
	if (loss_change_rate < cfg->em_converge_thresh) {
		if (cfg->ga_reset_sampling && sampling_method != cfg->ga_sampling_method2) {
			cfg->ga_reset_sampling = false;
			sampling_method        = cfg->ga_sampling_method2;
			std::cout << "Switched to sampling method: " + std::to_string(sampling_method) << std::endl;
		} else if (learning_rate > cfg->ending_step_size)
			learning_rate = std::max(learning_rate * 0.5f, cfg->ending_step_size);
		else
			count_no_progress += 1;
	} else
		count_no_progress = 0;
}

double EmModel::computeAndSyncLoss(std::vector<MolData> &data, suft_counts_t &suft, unsigned int energy) {

	double loss = 0.0;
	auto mol_it = data.begin();
	for (int molidx = 0; mol_it != data.end(); ++mol_it, molidx++) {
		if (mol_it->getGroup() != validation_group) {
			double mol_loss = computeLogLikelihoodLoss(molidx, *mol_it, suft, energy);
			loss += mol_loss;
		}
	}

	// update L2 only if lambda > 0
	if (cfg->lambda > 0.0) loss += getRegularizationTerm(energy);

	return loss;
}

double EmModel::gaGetUpdatedLearningRate(double learning_rate, int iter) const {

	if (USE_NO_DECAY == cfg->ga_decay_method) return learning_rate;

	if (USE_DEFAULT_DECAY == cfg->ga_decay_method)
		learning_rate *= 1.0 / (1.0 + cfg->decay_rate * (iter - 1));
	else if (USE_EXP_DECAY == cfg->ga_decay_method)
		learning_rate *= exp(-cfg->exp_decay_k * (iter - 1));
	else if (USE_STEP_DECAY == cfg->ga_decay_method)
		learning_rate *= pow(cfg->step_decay_drop, floor(iter / cfg->step_decay_epochs_drop));

	return learning_rate;
}

int EmModel::computeAndAccumulateGradient(double *grads, int mol_idx, MolData &mol_data, suft_counts_t &suft,
                                          int sampling_method, unsigned int energy) {

	unsigned int num_transitions = mol_data.getNumTransitions();
	unsigned int num_fragments   = mol_data.getNumFragments();

	int offset               = num_transitions;
	int num_used_transitions = 0;
	if (!mol_data.hasComputedGraph()) return num_used_transitions;

	suft_t *suft_values = &(suft.values[mol_idx]);

	unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();
	unsigned int suft_offset = energy * (num_transitions + num_fragments);

	std::set<int> selected_trans_id;
	if (sampling_method != USE_NO_SAMPLING) {
		getSubSampledTransitions(mol_data, sampling_method, energy, selected_trans_id);
		num_used_transitions = selected_trans_id.size();
	}

	// Iterate over from_id (i)
	auto frag_trans_map = mol_data.getFromIdTMap()->begin();
	for (int from_idx = 0; frag_trans_map != mol_data.getFromIdTMap()->end(); ++frag_trans_map, from_idx++) {

		// Do some random selection
		std::vector<int> sampled_ids;
		if (sampling_method != USE_NO_SAMPLING) {
			for (auto &id : *frag_trans_map)
				if (selected_trans_id.find(id) != selected_trans_id.end()) sampled_ids.push_back(id);
		} else
			sampled_ids = *frag_trans_map;

		// Calculate the denominator of the sum terms
		double denom = 1.0;
		for (auto trans_id : sampled_ids) denom += exp(mol_data.getThetaForIdx(energy, trans_id));

		// Complete the innermost sum terms	(sum over j')
		// std::map<unsigned int, double> sum_terms;
		std::vector<double> sum_terms(mol_data.getFeatureVectorForIdx(0)->getTotalLength());

		for (auto trans_id : sampled_ids) {
			const FeatureVector *fv = mol_data.getFeatureVectorForIdx(trans_id);

			for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
				auto fv_idx = *fv_it;
				double val  = exp(mol_data.getThetaForIdx(energy, trans_id)) / denom;
				sum_terms[fv_idx] += val;
			}
		}

		// Accumulate the transition (i \neq j) terms of the gradient (sum over j)
		double nu_sum = 0.0;
		for (auto trans_id : sampled_ids) {
			double nu = (*suft_values)[trans_id + suft_offset];
			nu_sum += nu;
			const FeatureVector *fv = mol_data.getFeatureVectorForIdx(trans_id);
			for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
				auto fv_idx = *fv_it;
				*(grads + fv_idx + grad_offset) += nu;
			}
		}

		// Accumulate the last term of each transition and the
		// persistence (i = j) terms of the gradient and Q
		double nu = (*suft_values)[offset + from_idx + suft_offset]; // persistence (i=j)
		for (auto idx = 0; idx < sum_terms.size(); ++idx)
			if (sum_terms[idx] != 0) *(grads + idx + grad_offset) -= (nu_sum + nu) * sum_terms[idx];
	}

	return num_used_transitions;
}

void EmModel::getSubSampledTransitions(MolData &moldata, int sampling_method, unsigned int energy,
                                       std::set<int> &selected_trans_id) const {

	switch (sampling_method) {
	case USE_RANDOM_SAMPLING: {
		moldata.getRandomSampledTransitions(selected_trans_id, cfg->ga_sampling_max_selection);
		break;
	}
	case USE_GRAPH_RANDOM_WALK_SAMPLING: {
		moldata.getSampledTransitionIdsRandomWalk(selected_trans_id, cfg->ga_sampling_max_selection);
		break;
	}
	case USE_DIFFERENCE_SAMPLING_BFS_CO:
	case USE_DIFFERENCE_SAMPLING_BFS: {
		// we are going to need  to predict spectrum without any postprocessing
		// apart from if we are using log scale
		moldata.computePredictedSpectra(*param, true, energy, 1, 1000, 100, 0.0, -1, cfg->use_log_scale_peak);

		std::set<unsigned int> selected_weights;

		moldata.getSelectedMasses(selected_weights, energy);
		if (sampling_method == USE_DIFFERENCE_SAMPLING_BFS_CO)
			moldata.getSampledTransitionIdUsingDiffMapBFS(selected_trans_id, selected_weights);
		else if (sampling_method == USE_DIFFERENCE_SAMPLING_BFS)
			moldata.getSampledTransitionIdUsingDiffMapCA(selected_trans_id, selected_weights);
		break;
	}
	default: break;
	}
}

double EmModel::computeLogLikelihoodLoss(int molidx, MolData &moldata, suft_counts_t &suft, unsigned int energy) {

	double q                     = 0.0;
	unsigned int num_transitions = moldata.getNumTransitions();
	unsigned int num_fragments   = moldata.getNumFragments();

	int offset = num_transitions;

	if (!moldata.hasComputedGraph()) return q;

	// Compute the latest transition thetas
	moldata.computeTransitionThetas(*param);
	suft_t *suft_values = &(suft.values[molidx]);

	// Compute
	unsigned int suft_offset = energy * (num_transitions + num_fragments);
	// Iterate over from_id (i)
	auto it                  = moldata.getFromIdTMap()->begin();
	for (int from_idx = 0; it != moldata.getFromIdTMap()->end(); ++it, from_idx++) {

		// Calculate the denominator of the sum terms
		double denom = 1.0;
		for (const auto &itt : *it) denom += exp(moldata.getThetaForIdx(energy, itt));

		// Accumulate the transition (i \neq j) terms of the gradient (sum over j)
		for (const auto &itt : *it) {
			double nu = (*suft_values)[itt + suft_offset];
			q += nu * (moldata.getThetaForIdx(energy, itt) - log(denom));
		}

		// Accumulate the last term of each transition and the
		// persistence (i = j) terms of the gradient and Q
		double nu = (*suft_values)[offset + from_idx + suft_offset]; // persistence (i=j)
		q -= nu * log(denom);
	}
	return q;
}

double EmModel::getRegularizationTerm(unsigned int energy) {

	double reg_term = 0.0;
	auto it         = this->used_idxs.begin();
	for (; it != this->used_idxs.end(); ++it) {
		double weight = param->getWeightAtIdx(*it);
		reg_term -= 0.5 * cfg->lambda * weight * weight;
	}

	// Remove the Bias terms (don't regularize the bias terms!)
	unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
	double bias                     = param->getWeightAtIdx(energy * weights_per_energy);
	reg_term += 0.5 * cfg->lambda * bias * bias;

	return reg_term;
}

void EmModel::updateGradientForRegularizationTerm(double *grads, unsigned int energy) {

	auto it = this->used_idxs.begin();
	for (; it != this->used_idxs.end(); ++it) {
		float weight = param->getWeightAtIdx(*it);
		*(grads + *it) -= cfg->lambda * weight;
	}

	// Remove the Bias terms (don't regularize the bias terms!)
	unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
	float bias                      = param->getWeightAtIdx(energy * weights_per_energy);
	*(grads + energy * weights_per_energy) += cfg->lambda * bias;
}
