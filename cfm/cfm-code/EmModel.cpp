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
#include "mpi.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include "EmModel.h"
#include "Comparators.h"

EmModel::EmModel(config_t *a_cfg, FeatureCalculator *an_fc,
                 std::string &a_status_filename, std::string initial_params_filename) {
    cfg = a_cfg;
    fc = an_fc;
    status_filename = a_status_filename;
    initComms();
    int num_energies_to_include = cfg->map_d_to_energy.back() + 1;
    if (initial_params_filename.empty()) {
        param = boost::shared_ptr<Param>(
                new Param(fc->getFeatureNames(), num_energies_to_include));
        initial_params_provided = false;
        comm->printToMasterOnly("EM: No initial params provided");
    } else {
        param = boost::shared_ptr<Param>(new Param(initial_params_filename));
        while (param->getNumEnergyLevels() < num_energies_to_include)
            param->appendRepeatedPrevEnergyParams();
        initial_params_provided = true;
        std::string msg =
                "EM: Initial params provided from " + initial_params_filename;
        comm->printToMasterOnly(msg.c_str());
    }
    sparse_params = true;
}

EmModel::~EmModel() { delete comm; }

void EmModel::computeThetas(MolData *moldata) {
    moldata->computeTransitionThetas(*param);
}

double
EmModel::trainModel(std::vector<MolData> &molDataSet, int group, std::string &out_param_filename, int energy_level) {

    unused_zeroed = 0;
    int iter = 0;

    if (!initial_params_provided)
        param->initWeights(cfg->param_init_type);
    comm->broadcastParamsOrigMpi(param.get());
    validation_group = group;

    // Write the initialised params to file (we may get want to reload and use
    // with saved suft state, even before updating)
    std::string init_out_param_filename = out_param_filename + "_init";
    if (comm->isMaster())
        writeParamsToFile(init_out_param_filename);

    // EM
    iter = 0;
    double loss;
    double prev_loss = -DBL_MAX;
    double best_loss = -DBL_MAX;

    // make of copy of learing rate
    // so we can share the save lr var over all em iterations
    //init some flags
    double learning_rate = cfg->starting_step_size;
    int sampling_method = cfg->ga_sampling_method;
    int em_no_progress_count = 0;
    bool switch_to_weighted_jaccard = false;

    molDataPreProcessing(molDataSet, energy_level);

    while (iter < MAX_EM_ITERATIONS) {
        std::string iter_out_param_filename =
                out_param_filename + "_" + std::to_string(iter);

        std::string msg = "EM Iteration " + std::to_string(iter);
        if (comm->isMaster())
            writeStatus(msg.c_str());
        comm->printToMasterOnly(msg.c_str());

        time_t before, after;
        std::vector<MolData>::iterator mol_it;

        // Reset sufficient counts
        suft_counts_t suft;
        initSuft(suft, molDataSet);

        int num_converged = 0, num_nonconverged = 0;
        int tot_numc = 0, total_numnonc = 0;
        before = time(nullptr);

        // Do the inference part (E-step)
        mol_it = molDataSet.begin();
        for (int molidx = 0; mol_it != molDataSet.end(); ++mol_it, molidx++) {

            if (!mol_it->hasComputedGraph())
                continue; // If we couldn't compute it's graph for some reason..
            if (mol_it->hasEmptySpectrum()) {
                std::cout << "Warning: No peaks with explanatory fragment found for "
                          << mol_it->getId() << ", ignoring this input molecule."
                          << std::endl;
                continue; // Ignore any molecule with poor (no peaks matched a fragment)
                // or missing spectra.
            }

            //if (mol_it->getGroup() == validation_group)
            //  continue;

            MolData *moldata = &(*mol_it);

            computeThetas(moldata);
            moldata->computeLogTransitionProbabilities();

            // Apply the peak evidence, compute the beliefs and record the sufficient
            // statistics

            if (cfg->use_single_energy_cfm) {
                beliefs_t beliefs;
                Inference infer(moldata, cfg);
                infer.calculateBeliefs(beliefs);
                recordSufficientStatistics(suft, molidx, moldata, &beliefs);
            } else {
                IPFP ipfp(moldata, cfg);
                beliefs_t *beliefs = ipfp.calculateBeliefs();
                int status = ipfp.status;
                if (status == NON_CONVERGE || status == OSC_CONVERGE)
                    num_nonconverged++;
                else if (status == COMPLETE_CONVERGE || status == CONVERGE_AFTER_MOD)
                    num_converged++;
                recordSufficientStatistics(suft, molidx, moldata, beliefs);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait
        after = time(nullptr);
        if (!cfg->use_single_energy_cfm) {
            total_numnonc = comm->collectSumInMaster(num_nonconverged);
            tot_numc = comm->collectSumInMaster(num_converged);
            std::string cvg_msg =
                    "Num Converged: " + std::to_string(tot_numc);
            std::string noncvg_msg = "Num Non-Converged: " +
                                     std::to_string(total_numnonc);
            if (comm->isMaster()) {
                writeStatus(cvg_msg.c_str());
                writeStatus(noncvg_msg.c_str());
            }
            comm->printToMasterOnly(cvg_msg.c_str());
            comm->printToMasterOnly(noncvg_msg.c_str());
        }
        std::string estep_time_msg =
                "[E-Step]Completed E-step processing: Time Elapsed = " +
                std::to_string(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(estep_time_msg.c_str());
        comm->printToMasterOnly(estep_time_msg.c_str());

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master
        // Find a new set of parameters to maximize the expected log likelihood
        // (M-step)

        before = time(nullptr);

        loss = updateParametersGradientAscent(molDataSet, suft, learning_rate, sampling_method);

        after = time(nullptr);
        std::string param_update_time_msg =
                "[M-Step]Completed M-step nn_param update: Time Elapsed = " +
                std::to_string(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(param_update_time_msg.c_str());
        comm->printToMasterOnly(param_update_time_msg.c_str());

        // Write the params
        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master

        before = time(nullptr);
        // validation Q value
        double val_q = 0.0;

        int molidx = 0, num_val_mols = 0, num_training_mols = 0;
        double jaccard = 0.0, w_jaccard = 0.0;
        for (mol_it = molDataSet.begin(); mol_it != molDataSet.end(); ++mol_it, molidx++) {
            if (mol_it->getGroup() == validation_group) {
                computeValidationMetrics(energy_level, molidx, mol_it, suft, val_q, num_val_mols, jaccard,
                                         w_jaccard);
            } else{
                /*double mol_loss = computeLogLikelihoodLoss(molidx, *mol_it, suft);
                std::cout << mol_it->getId() << " Loss=" << mol_loss;

                mol_it->computePredictedSpectra(*param, true, true, energy_level);
                Comparator *jaccard_cmp = new Jaccard(cfg->ppm_mass_tol, cfg->abs_mass_tol);
                Comparator *weighted_jaccard_cmp = new WeightedJaccard(cfg->ppm_mass_tol, cfg->abs_mass_tol);
                std::cout << " Jaccard="
                << jaccard_cmp->computeScore(mol_it->getSpectrum(energy_level),
                        mol_it->getPredictedSpectrum(energy_level))
                << " Weighted Jaccard=" << weighted_jaccard_cmp->computeScore(
                        mol_it->getSpectrum(energy_level),
                        mol_it->getPredictedSpectrum(energy_level)) << " ";
                delete jaccard_cmp;
                delete weighted_jaccard_cmp;
                std::cout << std::endl;*/
                num_training_mols++;
            }

        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master
        after = time(nullptr);
        std::string q_time_msg =
                "[M-Step] Finished Q compute: Time Elapsed = " +
                std::to_string(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(q_time_msg.c_str());

        val_q = comm->collectQInMaster(val_q);
        num_val_mols = comm->collectSumInMaster(num_val_mols);
        num_training_mols = comm->collectSumInMaster(num_training_mols);
        jaccard = comm->collectQInMaster(jaccard);
        w_jaccard = comm->collectQInMaster(w_jaccard);

        // Check for convergence
        double q_ratio = fabs((loss - prev_loss) / loss);
        double best_q_ratio = fabs((loss - best_loss) / loss);
        if (comm->isMaster()) {
            std::string qdif_str = "[M-Step]";
            if (prev_loss != -DBL_MAX)
                qdif_str += "Q_ratio= " + std::to_string(q_ratio) + " prev_Q=" + std::to_string(prev_loss) + "\n";

            if (best_loss != -DBL_MAX)
                qdif_str += "Best_Q_ratio= " + std::to_string(best_q_ratio) +
                            " best_Q=" + std::to_string(best_loss) + "\n";

            qdif_str += "Q=" + std::to_string(loss) + " Validation_Q=" + std::to_string(val_q) + "\n";
            qdif_str +=
            "Q_avg=" + std::to_string(loss / num_training_mols)
            + " Validation_Q_avg=" + std::to_string(val_q / num_val_mols)
            + " Validation Jaccard_Avg=" + std::to_string(jaccard / num_val_mols)
            + " Weighted Validation Jaccard_Avg=" += std::to_string(w_jaccard / num_val_mols);
            writeStatus(qdif_str.c_str());
            comm->printToMasterOnly(qdif_str.c_str());
        }

        // first let us save the model
        // only save the one has best Q so far
        if (best_loss < loss) {
            best_loss = loss;
            // Write the params
            if (comm->isMaster()) {
                std::string progress_str = "[M-Step] Found Better Q: "
                                           + std::to_string(best_loss) + " Write to File";
                comm->printToMasterOnly(progress_str.c_str());
                writeParamsToFile(iter_out_param_filename);
                writeParamsToFile(out_param_filename);
            }
        }

        // check if EM meet halt flag
        updateTraningParams(loss, prev_loss, q_ratio, learning_rate, sampling_method, em_no_progress_count);

        if (em_no_progress_count >= 3) {
            comm->printToMasterOnly(
                    ("EM Stopped after " + std::to_string(em_no_progress_count) + " No Progress Iterations").c_str());
            comm->printToMasterOnly(("EM Converged after " + std::to_string(iter) + " iterations").c_str());
            break;
        }

        prev_loss = loss;

        if (comm->isMaster()) {
            double loss_before_reg = loss - getRegularizationTerm();
            updateWJaccardFlag(switch_to_weighted_jaccard, prev_loss,
                               best_loss, loss_before_reg / num_training_mols, cfg->em_wjaccard_swicth_threshold);

        }
        switch_to_weighted_jaccard = comm->broadcastBooleanFlag(switch_to_weighted_jaccard);
        iter++;
    }

    if (iter >= MAX_EM_ITERATIONS)
        comm->printToMasterOnly(("Warning: EM did not converge after " +
                                 std::to_string(iter) +
                                 " iterations.")
                                        .c_str());

    return best_loss;
}

void EmModel::molDataPreProcessing(std::vector<MolData> &molDataSet, int energy_level) const {// pre process data
    for (auto &mol : molDataSet) {
        if (mol.getGroup() == validation_group)
            continue;

        mol.removePeaksWithNoFragment(cfg->abs_mass_tol, cfg->ppm_mass_tol);
        if (cfg->use_graph_pruning) {
            if (cfg->use_single_energy_cfm) {
                mol.createNewGraphForComputation();
                mol.pruneGraphBySpectra(energy_level, cfg->abs_mass_tol, cfg->ppm_mass_tol,
                                        cfg->aggressive_graph_pruning);
            } else {
                mol.pruneGraphBySpectra(-1, cfg->abs_mass_tol, cfg->ppm_mass_tol, cfg->aggressive_graph_pruning);
            }
        }
        if (cfg->add_noise) {
            mol.addNoise(cfg->noise_max, cfg->noise_sum, cfg->abs_mass_tol, cfg->ppm_mass_tol);
            mol.removePeaksWithNoFragment(cfg->abs_mass_tol, cfg->ppm_mass_tol);
        }
    }
}

void
EmModel::updateWJaccardFlag(bool &switch_to_weighted_jaccard, double &prev_loss, double &best_loss, double avg_loss,
                            double threshold) const {
    if (!switch_to_weighted_jaccard && avg_loss > threshold && cfg->em_use_weighted_jaccard) {
        switch_to_weighted_jaccard = true;
        prev_loss = 0.0;
        best_loss = 0.0;
        std::cout << "[EM INFO]Switching to Jaccard " << std::endl;
    }
}

void
EmModel::updateTraningParams(double loss, double prev_loss, double q_ratio, double &learning_rate, int &sampling_method,
                             int &count_no_progress) const {
    if (q_ratio < cfg->em_converge_thresh || prev_loss >= loss) {

        if (learning_rate > cfg->starting_step_size) {
            learning_rate *= 0.5;
            count_no_progress = 0;
        } else if (sampling_method != USE_NO_SAMPLING && cfg->reset_sampling) {
            if (comm->isMaster())
                std::cout << "[Reset] Turn off sampling" << std::endl;

            sampling_method = USE_NO_SAMPLING;
            learning_rate = cfg->starting_step_size * cfg->reset_sampling_lr_ratio;
            count_no_progress = 0;
        } else
            count_no_progress += 1;
    } else
        count_no_progress = 0;
}

double
EmModel::computeLoss(std::vector<MolData> &data, suft_counts_t &suft) {

    double loss = 0.0;
    auto mol_it = data.begin();
    for (int molidx = 0; mol_it != data.end(); ++mol_it, molidx++) {
        if (mol_it->getGroup() != validation_group) {
            double mol_loss = computeLogLikelihoodLoss(molidx, *mol_it, suft);
            loss += mol_loss;
        }
    }

    // update L2 only if lambda > 0
    if (comm->isMaster() && cfg->lambda > 0.0)
        loss += getRegularizationTerm();
    loss = comm->collectQInMaster(loss);
    loss = comm->broadcastQ(loss);

    return loss;
}

void EmModel::computeValidationMetrics(int energy_level, int molidx,
                                       std::vector<MolData, std::allocator<MolData>>::iterator &moldata,
                                       suft_counts_t &suft, double &val_q, int &num_val_mols, double &jaccard,
                                       double &w_jaccard) {

    val_q += computeLogLikelihoodLoss(molidx, *moldata, suft);
    num_val_mols++;
    Comparator *jaccard_cmp = new Jaccard(cfg->ppm_mass_tol, cfg->abs_mass_tol);
    Comparator *weighed_jaccard_cmp = new WeightedJaccard(cfg->ppm_mass_tol, cfg->abs_mass_tol);

    if (cfg->use_single_energy_cfm) {
        moldata->computePredictedSpectra(*param, false, false, energy_level);
        moldata->postprocessPredictedSpectra(100, 1, 30, 2.0);
        jaccard += jaccard_cmp->computeScore(moldata->getSpectrum(energy_level), moldata->getPredictedSpectrum(energy_level));
        w_jaccard += weighed_jaccard_cmp->computeScore(moldata->getSpectrum(energy_level),
                                          moldata->getPredictedSpectrum(energy_level));
    } else {
        moldata->computePredictedSpectra(*param, false, false);
        moldata->postprocessPredictedSpectra(100, 1, 30, 2.0);
        std::vector<unsigned int> energies;
        getEnergiesLevels(energies);
        for (auto &energy: energies) {
            jaccard += jaccard_cmp->computeScore(moldata->getSpectrum(energy), moldata->getPredictedSpectrum(energy));
            w_jaccard += weighed_jaccard_cmp->computeScore(moldata->getSpectrum(energy), moldata->getPredictedSpectrum(energy));
        }
        jaccard /= (double) energies.size();
        w_jaccard /= (double) energies.size();
    }
    delete jaccard_cmp;
    delete weighed_jaccard_cmp;
}

void EmModel::initSuft(suft_counts_t &suft, std::vector<MolData> &data) {

    // Resize the suft structure for each molecule
    auto num_mols = data.size();
    suft.values.resize(num_mols);
    for (unsigned int i = 0; i < num_mols; i++) {
        //const FragmentGraph *fg = data[i].getFragmentGraph();
        int len = data[i].getNumTransitions() + data[i].getNumFragments();
        int num_spectra = data[i].getNumSpectra();
        suft.values[i].resize(len * num_spectra);
    }
}

void EmModel::recordSufficientStatistics(suft_counts_t &suft, int molidx,
                                         MolData *moldata, beliefs_t *beliefs) {

    //const FragmentGraph *fg = moldata->getFragmentGraph();

    unsigned int num_transitions = moldata->getNumTransitions();
    unsigned int num_fragments = moldata->getNumFragments();

    int len_offset = num_transitions + num_fragments;

    // Accumulate the Sufficient Statistics
    for (unsigned int i = 0; i < num_transitions; i++) {

        const Transition *t = moldata->getTransitionAtIdx(i);

        double belief = 0.0;
        int energy = cfg->map_d_to_energy[0];
        if (t->getFromId() == 0) // main ion is always id = 0
            belief += exp(beliefs->tn[i][0]);

        for (unsigned int d = 1; d < cfg->model_depth; d++) {
            energy = cfg->map_d_to_energy[d];
            if (energy != cfg->map_d_to_energy[d - 1]) {
                suft.values[molidx][i + cfg->map_d_to_energy[d - 1] * len_offset] =
                        belief;
                belief = 0.0;
            }
            belief += exp(beliefs->tn[i][d]);
        }
        suft.values[molidx][i + energy * len_offset] = belief;
    }

    // Accumulate the persistence terms
    int offset = num_transitions;
    for (unsigned int i = 0; i < num_fragments; i++) {

        double belief = 0.0;
        int energy = cfg->map_d_to_energy[0];

        if (i == 0) // main ion is always id = 0
            belief += exp(beliefs->ps[i][0]);
        for (unsigned int d = 1; d < cfg->model_depth; d++) {
            energy = cfg->map_d_to_energy[d];
            if (energy != cfg->map_d_to_energy[d - 1]) {
                suft.values[molidx][i + offset +
                                    cfg->map_d_to_energy[d - 1] * len_offset] = belief;
                belief = 0.0;
            }
            belief += exp(beliefs->ps[i][d]);
        }
        suft.values[molidx][i + offset + energy * len_offset] = belief;
    }
}

double EmModel::updateParametersGradientAscent(std::vector<MolData> &data, suft_counts_t &suft, double learning_rate,
                                               int sampling_method) {

    // DBL_MIN is the smallest positive double
    // -DBL_MAX is the smallest negative double
    double loss = 0.0, prev_loss = -DBL_MAX, best_q = -DBL_MAX;

    std::vector<double> grads(param->getNumWeights(), 0.0);

    Solver *solver = nullptr;
    if (comm->isMaster()) {
        solver = getSolver(cfg->ga_method, learning_rate);
    }

    //Initial Q and gradient calculation (to determine used indexes)
    if (comm->used_idxs.empty()) {
        if (comm->isMaster())
            std::cout << "[M-Step] Initial Calculation ...";
        auto itdata = data.begin();
        for (int molidx = 0; itdata != data.end(); ++itdata, molidx++) {
            if (itdata->getGroup() != validation_group)
                computeAndAccumulateGradient(&grads[0], molidx, *itdata, suft, true, comm->used_idxs, 0);
        }

        comm->setMasterUsedIdxs();
        if (comm->isMaster())
            zeroUnusedParams();
        if (comm->isMaster())
            std::cout << "Done" << std::endl;
    }

    int n = 0;
    if (comm->isMaster())
        n = ((MasterComms *) comm)->master_used_idxs.size();
    comm->broadcastNumUsed(n);

    int iter = 0;
    //double learn_mult = 1.0;

    int max_iteration = cfg->ga_max_iterations;
    int ga_no_progress_count = 0;
    while (iter++ < max_iteration
           && fabs((loss - prev_loss) / loss) >= cfg->ga_converge_thresh
           && ga_no_progress_count <= 3) {

        if (iter > 1)
            prev_loss = loss;

        // update learning rate
        if (comm->isMaster() && USE_NO_DECAY != cfg->ga_decay_method) {
            learning_rate = getUpdatedLearningRate(learning_rate, loss, prev_loss, iter);
            solver->setLearningRate(learning_rate);
        }

        // Select molecules to include in gradient mini-batch.
        std::vector<int> minibatch_flags(data.size());
        int num_batch = cfg->ga_minibatch_nth_size;
        setMiniBatchFlags(minibatch_flags, num_batch);

        // Compute the gradient
        std::fill(grads.begin(), grads.end(), 0.0);
        auto mol_it = data.begin();
        for (auto batch_idx = 0; batch_idx < num_batch; ++batch_idx) {

            for (int molidx = 0; mol_it != data.end(); ++mol_it, molidx++) {
                if (minibatch_flags[molidx] == batch_idx && mol_it->getGroup() != validation_group) {
                    computeAndAccumulateGradient(&grads[0], molidx, *mol_it, suft, false, comm->used_idxs,
                                                 sampling_method);
                }
            }
            comm->collectGradsInMaster(&grads[0]);

            // Step the parameters
            if (comm->isMaster()) {
                // update L2 only if lambda > 0
                if(cfg->lambda > 0.0)
                    updateGradientForRegularizationTerm(&grads[0]);
                solver->adjustWeights(grads, ((MasterComms *) comm)->master_used_idxs, param);
            }

            // this should be a better way in large number of cores
            comm->broadcastParamsOrigMpi(param.get());
        }

        // End of batch
        // compute loss
        loss = computeLoss(data, suft);
        if (comm->isMaster()) {
            std::cout << iter << ":  Loss=" << loss << " Prev_Loss=" << prev_loss << " Learning_Rate=" << learning_rate
                      << std::endl;
            // let us roll Dropouts
            param->updateDropoutsRate(cfg->ga_dropout_delta, cfg->ga_dropout_lowerbond);
            param->rollDropouts(iter, cfg->ga_dropout_delta);
        }

        ga_no_progress_count = prev_loss >= loss ? ga_no_progress_count + 1 : 0;
    }

    if (comm->isMaster()) {
        if (iter == cfg->ga_max_iterations)
            std::cout << "Gradient ascent did not converge" << std::endl;
        else
            std::cout << "Gradient ascent converged after " << iter << " iterations"
                      << std::endl;
        delete solver;
    }
    return loss;
}

double EmModel::getUpdatedLearningRate(double learning_rate, double current_loss, double prev_loss, int iter) const {

    if (USE_NO_DECAY == cfg->ga_decay_method)
        return learning_rate;

    if (USE_DEFAULT_DECAY == cfg->ga_decay_method)
        learning_rate *= 1.0 / (1.0 + cfg->decay_rate * (iter - 1));
    else if (USE_EXP_DECAY == cfg->ga_decay_method)
        learning_rate *= exp(-cfg->exp_decay_k * iter);
    else if (USE_STEP_DECAY == cfg->ga_decay_method)
        learning_rate *= pow(cfg->step_decay_drop, floor(iter / cfg->step_decay_epochs_drop));

    if (current_loss < prev_loss && iter > 1 && learning_rate > cfg->starting_step_size * 0.02) {
        learning_rate = learning_rate * 0.5;
    }
    return learning_rate;
}

void EmModel::computeAndAccumulateGradient(double *grads, int mol_idx, MolData &mol_data, suft_counts_t &suft,
                                           bool record_used_idxs_only, std::set<unsigned int> &used_idxs,
                                           int sampling_method) {

    unsigned int num_transitions = mol_data.getNumTransitions();
    unsigned int num_fragments = mol_data.getNumFragments();

    int offset = num_transitions;

    if (!mol_data.hasComputedGraph())
        return;

    suft_t *suft_values = &(suft.values[mol_idx]);

    // Collect energies to compute
    std::vector<unsigned int> energies;
    getEnergiesLevels(energies);

    // Compute the gradients
    for (auto energy : energies) {
        unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();
        unsigned int suft_offset = energy * (num_transitions + num_fragments);

        std::set<int> selected_trans_id;
        if (!record_used_idxs_only && sampling_method != USE_NO_SAMPLING)
            getSubSampledTransitions(mol_data, sampling_method, energy, selected_trans_id);

        // Iterate over from_id (i)
        auto frag_trans_map = mol_data.getFromIdTMap()->begin();
        for (int from_idx = 0; frag_trans_map != mol_data.getFromIdTMap()->end(); ++frag_trans_map, from_idx++) {

            if (record_used_idxs_only) {
                for (auto trans_id : *frag_trans_map) {
                    const FeatureVector *fv = mol_data.getFeatureVectorForIdx(trans_id);
                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it)
                        used_idxs.insert(*fv_it + grad_offset);
                }
            } else {
                //Do some random selection
                std::vector<int> sampled_ids;
                if (sampling_method != USE_NO_SAMPLING){
                    for (auto id: *frag_trans_map)
                        if (selected_trans_id.find(id) != selected_trans_id.end())
                            sampled_ids.push_back(id);
                }
                else
                    sampled_ids = *frag_trans_map;

                // Calculate the denominator of the sum terms
                double denom = 1.0;
                for (auto trans_id : sampled_ids)
                    denom += exp(mol_data.getThetaForIdx(energy, trans_id));

                // Complete the innermost sum terms	(sum over j')
                std::map<unsigned int, double> sum_terms;

                for (auto trans_id : sampled_ids) {
                    const FeatureVector *fv = mol_data.getFeatureVectorForIdx(trans_id);

                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                        auto fv_idx = *fv_it;
                        double val = exp(mol_data.getThetaForIdx(energy, trans_id)) / denom;
                        if (sum_terms.find(fv_idx) != sum_terms.end())
                            sum_terms[fv_idx] += val;
                        else
                            sum_terms[fv_idx] = val;
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
                for (auto &sit: sum_terms) {
                    *(grads + sit.first + grad_offset) -= (nu_sum + nu) * sit.second;
                }
            }
        }
    }

    // Compute the latest transition thetas
    if (!record_used_idxs_only)
        mol_data.computeTransitionThetas(*param);
}


void EmModel::getSubSampledTransitions(MolData &moldata, int sampling_method, unsigned int energy,
                                       std::set<int> &selected_trans_id) const {

    int num_trans = moldata.getNumTransitions();
    int num_iterations = (int) ((cfg->ga_graph_sampling_k * num_trans) / (double) (cfg->fg_depth * cfg->fg_depth));
    if (sampling_method == USE_GRAPH_WEIGHTED_RANDOM_WALK_SAMPLING) {
        moldata.computePredictedSpectra(*param, false, false, energy);
        moldata.getSampledTransitionIdsWeightedRandomWalk(selected_trans_id, num_iterations, energy,
                                                          moldata.getWeightedJaccardScore(energy));
    } else if (sampling_method == USE_GRAPH_RANDOM_WALK_SAMPLING) {
        moldata.getSampledTransitionIdsRandomWalk(selected_trans_id, 0.1);
    } else if (sampling_method == USE_DIFFERENCE_SAMPLING) {
        moldata.computePredictedSpectra(*param, false, false, energy);
        std::set<double> selected_weights;
        std::set<double> all_weights;

        moldata.getSelectedWeights(selected_weights, all_weights, energy, true);
        moldata.getSampledTransitionIdUsingDiffMap(selected_trans_id, all_weights, all_weights);
        // std::cout << selected_trans_id.size() << "/" << moldata.getNumTransitions() << std::endl;
    }
}

double EmModel::computeLogLikelihoodLoss(int molidx, MolData &moldata, suft_counts_t &suft) {

    double q = 0.0;
    unsigned int num_transitions = moldata.getNumTransitions();
    unsigned int num_fragments = moldata.getNumFragments();

    int offset = num_transitions;

    if (!moldata.hasComputedGraph())
        return q;

    // Compute the latest transition thetas
    moldata.computeTransitionThetas(*param);
    suft_t *suft_values = &(suft.values[molidx]);

    // Collect energies to compute
    std::vector<unsigned int> energies;
    getEnergiesLevels(energies);

    // Compute
    for (auto energy : energies) {
        unsigned int suft_offset = energy * (num_transitions + num_fragments);

        // Iterate over from_id (i)
        auto it = moldata.getFromIdTMap()->begin();
        for (int from_idx = 0; it != moldata.getFromIdTMap()->end(); ++it, from_idx++) {

            // Calculate the denominator of the sum terms
            double denom = 1.0;
            for (auto itt : *it)
                denom += exp(moldata.getThetaForIdx(energy, itt));

            // Accumulate the transition (i \neq j) terms of the gradient (sum over j)
            for (auto itt : *it) {
                double nu = (*suft_values)[itt + suft_offset];
                q += nu * (moldata.getThetaForIdx(energy, itt) - log(denom));
            }

            // Accumulate the last term of each transition and the
            // persistence (i = j) terms of the gradient and Q
            double nu = (*suft_values)[offset + from_idx + suft_offset]; // persistence (i=j)
            q -= nu * log(denom);
        }
    }

    return q;
}


double EmModel::getRegularizationTerm() {

    double reg_term = 0.0;
    auto it = ((MasterComms *) comm)->master_used_idxs.begin();
    for (; it != ((MasterComms *) comm)->master_used_idxs.end(); ++it) {
        double weight = param->getWeightAtIdx(*it);
        reg_term -= 0.5 * cfg->lambda * weight * weight;
    }

    // Remove the Bias terms (don't regularize the bias terms!)
    unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
    for (unsigned int energy = 0; energy < param->getNumEnergyLevels();
         energy++) {
        double bias = param->getWeightAtIdx(energy * weights_per_energy);
        reg_term += 0.5 * cfg->lambda * bias * bias;
    }
    return reg_term;
}

void EmModel::updateGradientForRegularizationTerm(double *grads) {

    auto it = ((MasterComms *) comm)->master_used_idxs.begin();
    for (; it != ((MasterComms *) comm)->master_used_idxs.end(); ++it) {
        double weight = param->getWeightAtIdx(*it);
        *(grads + *it) -= cfg->lambda * weight;
    }

    // Remove the Bias terms (don't regularize the bias terms!)
    unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
    for (unsigned int energy = 0; energy < param->getNumEnergyLevels();
         energy++) {
        double bias = param->getWeightAtIdx(energy * weights_per_energy);
        *(grads + energy * weights_per_energy) += cfg->lambda * bias;
    }
}