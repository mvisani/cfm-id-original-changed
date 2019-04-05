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

    // mpi broadcast init weights and dropouts for nn
    comm->broadcastParamsWeightsOrigMpi(param.get());
    comm->broadcastDropouts(param.get());

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
    float learning_rate = cfg->starting_step_size;
    int sampling_method = cfg->ga_sampling_method;
    int em_no_progress_count = 0;

    auto mol_it = molDataSet.begin();
    for (int molidx = 0; mol_it != molDataSet.end(); ++mol_it, molidx++) {
        if (mol_it->hasEmptySpectrum(energy_level) && mol_it->getGroup() != validation_group)
            std::cout << "Warning: No peaks with explanatory fragment found for "
                      << mol_it->getId() << ", ignoring this input molecule."
                      << std::endl;
    }

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
            if (mol_it->hasEmptySpectrum(energy_level))
                continue; // Ignore any molecule with poor (no peaks matched a fragment)
            // or missing spectra.

            // do noting if disbaled
            if (mol_it->getGroup() == validation_group && cfg->disable_cross_val_metrics)
                continue;

            MolData *moldata = &(*mol_it);

            computeThetas(moldata);
            moldata->computeLogTransitionProbabilities();

            // Apply the peak evidence, compute the beliefs and record the sufficient
            // statistics

            beliefs_t beliefs;
            Inference infer(moldata, cfg);
            infer.calculateBeliefs(beliefs, energy_level);
            recordSufficientStatistics(suft, molidx, moldata, &beliefs, energy_level);
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait
        after = time(nullptr);

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

        loss = updateParametersGradientAscent(molDataSet, suft, learning_rate, sampling_method, energy_level);

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
        double train_jaccard = 0.0, train_w_jaccard = 0.0;
        double val_jaccard = 0.0, val_w_jaccard = 0.0;
        for (mol_it = molDataSet.begin(); mol_it != molDataSet.end(); ++mol_it, molidx++) {
            if (mol_it->getGroup() == validation_group && !cfg->disable_cross_val_metrics){
                num_val_mols++;
                val_q += computeLogLikelihoodLoss(molidx, *mol_it, suft, energy_level);
                computeMetrics(energy_level, mol_it, val_jaccard,
                               val_w_jaccard);

            }
            else{
                num_training_mols++;
                if(!cfg->disable_training_metrics){
                    computeMetrics(energy_level, mol_it, train_jaccard,
                                   train_w_jaccard);
                }
            }

        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master
        after = time(nullptr);
        std::string q_time_msg =
                "[M-Step] Finished Q compute: Time Elapsed = " +
                std::to_string(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(q_time_msg.c_str());

        num_training_mols = comm->collectSumInMaster(num_training_mols);

        if (!cfg->disable_training_metrics) {
            train_jaccard = comm->collectQInMaster(train_jaccard);
            train_w_jaccard = comm->collectQInMaster(train_w_jaccard);
        }

        if (!cfg->disable_cross_val_metrics) {
            val_q = comm->collectQInMaster(val_q);
            num_val_mols = comm->collectSumInMaster(num_val_mols);
            val_jaccard = comm->collectQInMaster(val_jaccard);
            val_w_jaccard = comm->collectQInMaster(val_w_jaccard);
        }

        // Check for convergence
        double loss_ratio = fabs((loss - prev_loss) / loss);
        if (comm->isMaster()) {
            std::string qdif_str = "[M-Step]";
            qdif_str += "Loss=" + std::to_string(loss) + " Loss_Avg=" + std::to_string(loss / num_training_mols);

            if (prev_loss != -DBL_MAX)
                qdif_str += "\nLoss_Ratio= " + std::to_string(loss_ratio) + " Prev_Loss=" + std::to_string(prev_loss);

            if (best_loss != -DBL_MAX)
                qdif_str += " Best_Loss=" + std::to_string(best_loss);

            if (!cfg->disable_training_metrics) {
                qdif_str += "\nJaccard_Avg=" + std::to_string(train_jaccard / num_training_mols)
                            + " Weighted_Jaccard_Avg=" += std::to_string(train_w_jaccard / num_training_mols);
            }

            if (!cfg->disable_cross_val_metrics) {
                qdif_str += "\nValidation_Loss=" + std::to_string(val_q)
                            + " Validation_Loss_Avg=" + std::to_string(val_q / num_val_mols)
                            + "\nValidation_Jaccard_Avg=" + std::to_string(val_jaccard / num_val_mols)
                            + " Weighted_Validation_Jaccard_Avg=" += std::to_string(val_w_jaccard / num_val_mols);
            }

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
        updateTrainingParams(loss, prev_loss, loss_ratio, learning_rate, sampling_method, em_no_progress_count);

        if (em_no_progress_count >= cfg->em_no_progress_count) {
            if (cfg->ga_reset_sampling && cfg->ga_sampling_method != cfg->ga_sampling_method2){
                cfg->ga_reset_sampling = false;
                cfg->ga_sampling_method = cfg->ga_sampling_method2;
            }
            else{
                comm->printToMasterOnly(
                        ("EM Stopped after " + std::to_string(em_no_progress_count) + " No Progress Iterations").c_str());
                comm->printToMasterOnly(("EM Converged after " + std::to_string(iter) + " iterations").c_str());
                //time to stop
                //before stop, let us load best model
                param->readFromFile(out_param_filename);
                break;
            }
        }

        prev_loss = loss;
        iter++;
    }

    if (iter >= MAX_EM_ITERATIONS)
        comm->printToMasterOnly(("Warning: EM did not converge after " +
                                 std::to_string(iter) +
                                 " iterations.")
                                        .c_str());

    return best_loss;
}

void
EmModel::updateTrainingParams(double loss, double prev_loss, double loss_ratio, float &learning_rate,
                              int &sampling_method,
                              int &count_no_progress) const {
    if (loss_ratio < cfg->em_converge_thresh || prev_loss >= loss) {
        if (learning_rate > cfg->ending_step_size)
            learning_rate = std::max(learning_rate * 0.5f, cfg->ending_step_size);
        else
            count_no_progress += 1;
    } else
        count_no_progress = 0;
}

double
EmModel::computeLoss(std::vector<MolData> &data, suft_counts_t &suft, unsigned int energy) {

    double loss = 0.0;
    auto mol_it = data.begin();
    for (int molidx = 0; mol_it != data.end(); ++mol_it, molidx++) {
        if (mol_it->getGroup() != validation_group) {
            double mol_loss = computeLogLikelihoodLoss(molidx, *mol_it, suft, energy);
            loss += mol_loss;
        }
    }

    // update L2 only if lambda > 0
    if (comm->isMaster() && cfg->lambda > 0.0)
        loss += getRegularizationTerm(energy);
    loss = comm->collectQInMaster(loss);
    loss = comm->broadcastQ(loss);

    return loss;
}

void EmModel::computeMetrics(int energy_level, std::vector<MolData, std::allocator<MolData>>::iterator &moldata,
                             double &jaccard, double &w_jaccard) {

    Comparator *jaccard_cmp = new Jaccard(cfg->ppm_mass_tol, cfg->abs_mass_tol);
    Comparator *weighed_jaccard_cmp = new WeightedJaccard(cfg->ppm_mass_tol, cfg->abs_mass_tol);

    moldata->computePredictedSpectra(*param, false, false, energy_level);
    moldata->postprocessPredictedSpectra(80, 1, 1000);
    jaccard += jaccard_cmp->computeScore(moldata->getOrigSpectrum(energy_level),
                                             moldata->getPredictedSpectrum(energy_level));
    w_jaccard += weighed_jaccard_cmp->computeScore(moldata->getOrigSpectrum(energy_level),
                                                       moldata->getPredictedSpectrum(energy_level));

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

void EmModel::recordSufficientStatistics(suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs,
                                         unsigned int energy) {

    int depth = cfg->model_depth > moldata->getFGHeight() ? cfg->model_depth : moldata->getFGHeight();
    unsigned int num_transitions = moldata->getNumTransitions();
    unsigned int num_fragments = moldata->getNumFragments();
    int len_offset = num_transitions + num_fragments;

    // Accumulate the Sufficient Statistics
    for (unsigned int i = 0; i < num_transitions; i++) {

        const TransitionPtr t = moldata->getTransitionAtIdx(i);

        double belief = 0.0;
        //int energy = cfg->map_d_to_energy[0];
        if (t->getFromId() == 0) // main ion is always id = 0
            belief += exp(beliefs->tn[i][0]);

        for (unsigned int d = 1; d < depth; d++) {
            /*energy = cfg->map_d_to_energy[d];
            if (energy != cfg->map_d_to_energy[d - 1]) {
                suft.values[molidx][i + cfg->map_d_to_energy[d - 1] * len_offset] =
                        belief;
                belief = 0.0;
            }*/
            belief += exp(beliefs->tn[i][d]);
        }
        suft.values[molidx][i + energy * len_offset] = belief;
    }

    // Accumulate the persistence terms
    int offset = num_transitions;
    for (unsigned int i = 0; i < num_fragments; i++) {

        double belief = 0.0;
        //int energy = cfg->map_d_to_energy[0];

        if (i == 0) // main ion is always id = 0
            belief += exp(beliefs->ps[i][0]);
        for (unsigned int d = 1; d < depth; d++) {
            /*energy = cfg->map_d_to_energy[d];
            if (energy != cfg->map_d_to_energy[d - 1]) {
                suft.values[molidx][i + offset +
                                    cfg->map_d_to_energy[d - 1] * len_offset] = belief;
                belief = 0.0;
            }*/
            belief += exp(beliefs->ps[i][d]);
        }
        suft.values[molidx][i + offset + energy * len_offset] = belief;
    }
}

double EmModel::updateParametersGradientAscent(std::vector<MolData> &data, suft_counts_t &suft, double learning_rate,
                                               int sampling_method, unsigned int energy) {

    // DBL_MIN is the smallest positive double
    // -DBL_MAX is the smallest negative double
    double loss = 0.0, prev_loss = -DBL_MAX, Best_Loss = -DBL_MAX;

    std::vector<float> grads(param->getNumWeights(), 0.0);

    Solver *solver = nullptr;
    if (comm->isMaster()) {
        solver = getSolver(cfg->ga_method, learning_rate);
    }

    //Initial Q and gradient calculation (to determine used indexes)
    if (comm->used_idxs.empty()) {
        if (comm->isMaster())
            std::cout << "[M-Step] Collect Used Index ...";

        auto itdata = data.begin();
        for (int molidx = 0; itdata != data.end(); ++itdata, molidx++) {
            if (itdata->getGroup() != validation_group){
                collectUsedIdx(*itdata, comm->used_idxs, energy);
                computeThetas(&(*itdata));
            }
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
    while (iter++ < max_iteration && ga_no_progress_count <= cfg->ga_no_progress_count) {

        if (iter > 1)
            prev_loss = loss;

        // update learning rate
        if (comm->isMaster() && USE_NO_DECAY != cfg->ga_decay_method) {
            learning_rate = getUpdatedLearningRate(learning_rate, iter);
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
            double num_trans = 0;
            for (int molidx = 0; mol_it != data.end(); ++mol_it, molidx++) {
                if (minibatch_flags[molidx] == batch_idx && mol_it->getGroup() != validation_group) {
                    // so now it should not crash anymore
                    if (!mol_it->hasComputedGraph())
                        continue;

                    num_trans += computeAndAccumulateGradient(&grads[0], molidx, *mol_it, suft,
                                                              sampling_method, energy);
                }
            }

            if(num_trans >0)
                for(auto & grad : grads)
                    grad /= num_trans;

            comm->collectGradsInMasterOrigMpi(grads);

            // Step the parameters
            if (comm->isMaster()) {
                // update L2 only if lambda > 0
                if (cfg->lambda > 0.0)
                    updateGradientForRegularizationTerm(&grads[0], energy);
                solver->adjustWeights(grads, ((MasterComms *) comm)->master_used_idxs, param);
            }

            // this should be a better way in large number of cores
            // comm->broadcastParamsWeights(param.get());
            comm->broadcastParamsWeightsOrigMpi(param.get());
        }

        // End of batch
        // compute loss
        loss = computeLoss(data, suft, energy);
        if (comm->isMaster()) {
            std::cout << iter << ":  Loss=" << loss << " Prev_Loss=" << prev_loss << " Learning_Rate=" << learning_rate
                      << std::endl;
            // let us roll Dropouts
            param->rollDropouts();
        }
        comm->broadcastDropouts(param.get());

        if((fabs((loss - prev_loss) / loss) < cfg->ga_converge_thresh) || (prev_loss >= loss))
            ga_no_progress_count++;
        else
            ga_no_progress_count = 0;
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

double EmModel::getUpdatedLearningRate(double learning_rate, int iter) const {

    if (USE_NO_DECAY == cfg->ga_decay_method)
        return learning_rate;

    if (USE_DEFAULT_DECAY == cfg->ga_decay_method)
        learning_rate *= 1.0 / (1.0 + cfg->decay_rate * (iter - 1));
    else if (USE_EXP_DECAY == cfg->ga_decay_method)
        learning_rate *= exp(-cfg->exp_decay_k * iter);
    else if (USE_STEP_DECAY == cfg->ga_decay_method)
        learning_rate *= pow(cfg->step_decay_drop, floor(iter / cfg->step_decay_epochs_drop));

    return learning_rate;
}

int EmModel::computeAndAccumulateGradient(float *grads, int mol_idx, MolData &mol_data, suft_counts_t &suft,
                                          int sampling_method, unsigned int energy) {

    unsigned int num_transitions = mol_data.getNumTransitions();
    unsigned int num_fragments = mol_data.getNumFragments();

    int offset = num_transitions;
    int num_used_transitions = 0;
    if (!mol_data.hasComputedGraph())
        return num_used_transitions;

    suft_t *suft_values = &(suft.values[mol_idx]);

    unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();
    unsigned int suft_offset = energy * (num_transitions + num_fragments);

    std::set<int> selected_trans_id;
    if (sampling_method != USE_NO_SAMPLING){
        getSubSampledTransitions(mol_data, sampling_method, energy, selected_trans_id);
        num_used_transitions = selected_trans_id.size();
    }

    // Iterate over from_id (i)
    auto frag_trans_map = mol_data.getFromIdTMap()->begin();
    for (int from_idx = 0; frag_trans_map != mol_data.getFromIdTMap()->end(); ++frag_trans_map, from_idx++) {

            //Do some random selection
            std::vector<int> sampled_ids;
            if (sampling_method != USE_NO_SAMPLING) {
                for (auto & id: *frag_trans_map)
                    if (selected_trans_id.find(id) != selected_trans_id.end())
                        sampled_ids.push_back(id);
            } else
                sampled_ids = *frag_trans_map;

            // Calculate the denominator of the sum terms
            double denom = 1.0;
            for (auto trans_id : sampled_ids)
                denom += exp(mol_data.getThetaForIdx(energy, trans_id));

            // Complete the innermost sum terms	(sum over j')
            //std::map<unsigned int, double> sum_terms;
            std::vector<double> sum_terms(mol_data.getFeatureVectorForIdx(0)->getTotalLength());

            for (auto trans_id : sampled_ids) {
                const FeatureVector *fv = mol_data.getFeatureVectorForIdx(trans_id);

                for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                    auto fv_idx = *fv_it;
                    double val = exp(mol_data.getThetaForIdx(energy, trans_id)) / denom;
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
            for(auto idx = 0 ; idx < sum_terms.size(); ++idx)
                if(sum_terms[idx] != 0)
                    *(grads + idx + grad_offset) -= (nu_sum + nu) * sum_terms[idx];
        }

    return num_used_transitions;
}

void EmModel::collectUsedIdx(MolData &mol_data, std::set<unsigned int> &used_idxs, unsigned int energy){

    if (!mol_data.hasComputedGraph())
        return;

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

void EmModel::getSubSampledTransitions(MolData &moldata, int sampling_method, unsigned int energy,
                                       std::set<int> &selected_trans_id) const {

    switch (sampling_method) {
        case USE_RANDOM_SAMPLING:
        {
            moldata.getRandomSampledTransitions(selected_trans_id, cfg->ga_sampling_max_selection);
            break;
        }
        case USE_GRAPH_RANDOM_WALK_SAMPLING:
        {
            moldata.getSampledTransitionIdsRandomWalk(selected_trans_id,cfg->ga_sampling_max_selection);
            break;
        }
        case USE_DIFFERENCE_SAMPLING:{
            moldata.computePredictedSpectra(*param, false, true, energy);
            std::set<unsigned int> selected_weights;

            moldata.getSelectedWeights(selected_weights, energy);
            moldata.getSampledTransitionIdUsingDiffMap(selected_trans_id, selected_weights);
            break;
        }
        default:
            break;
    }
}

double EmModel::computeLogLikelihoodLoss(int molidx, MolData &moldata, suft_counts_t &suft, unsigned int energy) {

    double q = 0.0;
    unsigned int num_transitions = moldata.getNumTransitions();
    unsigned int num_fragments = moldata.getNumFragments();

    int offset = num_transitions;

    if (!moldata.hasComputedGraph())
        return q;

    // Compute the latest transition thetas
    moldata.computeTransitionThetas(*param);
    suft_t *suft_values = &(suft.values[molidx]);

    // Compute
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
    return q;
}


double EmModel::getRegularizationTerm(unsigned int energy) {

    double reg_term = 0.0;
    auto it = ((MasterComms *) comm)->master_used_idxs.begin();
    for (; it != ((MasterComms *) comm)->master_used_idxs.end(); ++it) {
        if(withinGradOffset(*it, energy)){
            double weight = param->getWeightAtIdx(*it);
            reg_term -= 0.5 * cfg->lambda * weight * weight;
        }
    }

    // Remove the Bias terms (don't regularize the bias terms!)
    unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
    for (unsigned int energy = 0; energy < param->getNumEnergyLevels();
         energy++) {
        if(withinGradOffset(*it, energy)) {
            double bias = param->getWeightAtIdx(energy * weights_per_energy);
            reg_term += 0.5 * cfg->lambda * bias * bias;
        }
    }
    return reg_term;
}

void EmModel::updateGradientForRegularizationTerm(float *grads, unsigned int energy) {

    auto it = ((MasterComms *) comm)->master_used_idxs.begin();
    for (; it != ((MasterComms *) comm)->master_used_idxs.end(); ++it) {
        if(withinGradOffset(*it, energy)) {
            float weight = param->getWeightAtIdx(*it);
            *(grads + *it) -= cfg->lambda * weight;
        }
    }

    // Remove the Bias terms (don't regularize the bias terms!)
    unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
    for (unsigned int energy = 0; energy < param->getNumEnergyLevels();
         energy++) {
        if(withinGradOffset(*it, energy)) {
            float bias = param->getWeightAtIdx(energy * weights_per_energy);
            *(grads + energy * weights_per_energy) += cfg->lambda * bias;
        }
    }
}