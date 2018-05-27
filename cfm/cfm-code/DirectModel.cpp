//
// Created by feiw on 26/05/18.
//

#include "DirectModel.h"

double DirectModel::ComputeAndAccumulateGradient(double *grads, MolData &moldata,
                                                 bool record_used_idxs_only, std::set<unsigned int> &used_idxs,
                                                 int sampling_method) {
    double Q = 0.0;
    const FragmentGraph *fg = moldata.GetFragmentGraph();
    unsigned int num_transitions = fg->getNumTransitions();
    unsigned int num_fragments = fg->getNumFragments();

    int offset = num_transitions;

    if (!moldata.HasComputedGraph())
        return Q;

    // Collect energies to compute
    std::vector<unsigned int> energies;
    getEnergiesLevels(energies);

    // Compute the gradients
    for (auto energy : energies) {
        unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();
        unsigned int suft_offset = energy * (num_transitions + num_fragments);

        std::set<int> selected_trans_id;
        if (!record_used_idxs_only && sampling_method == USE_GRAPH_RANDOM_WALK_SAMPLING) {

            int num_frags = moldata.GetSpectrum(energy)->size();
            int num_iterations = cfg->ga_graph_sampling_k * num_frags;
            if (cfg->ga_use_sqaured_iter_num)
                num_iterations = num_iterations * (energy + 1) * (energy + 1);
            moldata.GetSampledTransitionIdsRandomWalk(selected_trans_id, num_iterations, energy, 1.0);
        }

        // Iterate over from_id (i)
        auto frag_trans_map = fg->getFromIdTMap()->begin();
        for (int from_idx = 0; frag_trans_map != fg->getFromIdTMap()->end(); ++frag_trans_map, from_idx++) {

            if (record_used_idxs_only) {
                for (auto trans_id : *frag_trans_map) {
                    const FeatureVector *fv = moldata.GetFeatureVectorForIdx(trans_id);
                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                        used_idxs.insert(*fv_it + grad_offset);
                    }
                }
            } else {
                //Do some random selection
                /*std::vector<int> sampled_ids;
                if (sampling_method == USE_GRAPH_RANDOM_WALK_SAMPLING) {
                    for (auto id: *frag_trans_map) {
                        if (selected_trans_id.find(id) != selected_trans_id.end()) {
                            sampled_ids.push_back(id);
                        }
                    }
                } else
                    sampled_ids = *frag_trans_map;



                // Calculate the denominator of the sum terms
                double denom = 1.0;
                for (auto trans_id : sampled_ids)
                    denom += exp(moldata.getThetaForIdx(energy, trans_id));

                // Complete the innermost sum terms	(sum over j')
                std::map<unsigned int, double> sum_terms;

                for (auto trans_id : sampled_ids) {
                    const FeatureVector *fv = moldata.getFeatureVectorForIdx(trans_id);

                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                        auto fv_idx = *fv_it;
                        double val = exp(moldata.getThetaForIdx(energy, trans_id)) / denom;
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
                    const FeatureVector *fv = moldata.GetFeatureVectorForIdx(trans_id);

                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                        auto fv_idx = *fv_it;
                        *(grads + fv_idx + grad_offset) += nu;

                        if (record_used_idxs_only)
                            used_idxs.insert(fv_idx + grad_offset);

                    }
                    Q += nu * (moldata.GetThetaForIdx(energy, trans_id) - log(denom));
                }

                // Accumulate the last term of each transition and the
                // persistence (i = j) terms of the gradient and Q
                double nu = (*suft_values)[offset + from_idx + suft_offset]; // persistence (i=j)
                for (auto &sit: sum_terms) {
                    *(grads + sit.first + grad_offset) -= (nu_sum + nu) * sit.second;

                    if (record_used_idxs_only)
                        used_idxs.insert(sit.first + grad_offset);
                }
                Q -= nu * log(denom);*/
            }
        }
    }

    // Compute the latest transition thetas
    if (!record_used_idxs_only)
        moldata.ComputeNormalizedTransitionThetas(*param);

    return Q;
}

double DirectModel::ComputeQ(MolData &moldata) {
    double Q = 0.0;
    moldata.ComputeLogTransitionProbabilities();
    // TODO FIX eng level
    int energy_level = 0;
    for(const auto & peak : *moldata.GetSpectrum(energy_level)){
        // get intensity or height of the specturm
        double intensity = peak.intensity;
        double path_log_prob = 0;
        // get possible path to this peak
        // anything out of mass tol range will not be considered
        std::vector<Path> pathes;
        double mass_tol = getMassTol(cfg->abs_mass_tol,cfg->ppm_mass_tol,peak.mass);
        moldata.GetPathes(pathes,peak.mass, mass_tol);

        int num_transitions = moldata.GetFragmentGraph()->getNumTransitions();
        for(auto & path: pathes){
            for(auto & trans_id : *path.GetTransIds()){
                if(trans_id < num_transitions){
                    Q += moldata.GetLogTransitionProbForIdx(energy_level, trans_id);
                } else {
                    Q += moldata.GetLogPersistenceProbForIdx(energy_level, trans_id);
                }
            }
        }
        Q *= intensity;
    }
    return Q;
}

double DirectModel::addRegularizersAndUpdateGradient(double *grads) {
    return 0.0;
}

double DirectModel::updateParametersGradientAscent(std::vector<MolData> &data,
                                                   double learning_rate, int sampling_method) {
    return 0.0;
}

double DirectModel::TrainModel(std::vector<MolData> &molDataSet, int group, std::string &out_param_filename) {
    unused_zeroed = 0;
    int iter = 0;

    if (!initial_params_provided)
        initParams();
    comm->broadcastInitialParams(param.get());
    validation_group = group;

    // Write the initialised params to file (we may get want to reload and use
    // with saved suft state, even before updating)
    std::string init_out_param_filename = out_param_filename + "_init";
    if (comm->isMaster())
        writeParamsToFile(init_out_param_filename);

    // Model
    iter = 0;
    double Q;
    double prevQ = -DBL_MAX;
    double bestQ = -DBL_MAX;

    // make of copy of learing rate
    // so we can share the save lr var over all em iterations
    double learning_rate = cfg->starting_step_size;
    int sampling_method = cfg->ga_sampling_method;

    int count_no_progress = 0;

    if (cfg->add_noise) {
        for (auto &mol : molDataSet) {
            if (mol.GetGroup() != validation_group)
                mol.AddNoise(cfg->noise_max, cfg->noise_sum, cfg->abs_mass_tol, cfg->ppm_mass_tol);
            mol.removePeaksWithNoFragment(cfg->abs_mass_tol, cfg->ppm_mass_tol);
        }
    }
    while (iter < 1000) {
        time_t before, after;
        std::string iter_out_param_filename =
                out_param_filename + "_" + boost::lexical_cast<std::string>(iter);

        std::string msg = "DirectModel Iteration " + boost::lexical_cast<std::string>(iter);

        int num_converged = 0, num_nonconverged = 0;
        int tot_numc = 0, total_numnonc = 0;

        Q = updateParametersGradientAscent(molDataSet, learning_rate, sampling_method);

        after = time(nullptr);
        std::string param_update_time_msg =
                "[M-Step]Completed M-step param update: Time Elapsed = " +
                boost::lexical_cast<std::string>(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(param_update_time_msg.c_str());
        comm->printToMasterOnly(param_update_time_msg.c_str());

        // Write the params
        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master

        before = time(nullptr);
        // validation Q value
        double valQ = 0.0;
        // Compute the final Q (with all molecules, in case only some were used in
        // the mini-batch)
        if (cfg->ga_minibatch_nth_size > 1 ||
            sampling_method != USE_NO_SAMPLING) {
            Q = 0;
            if (comm->isMaster())
                std::cout << "[INFO]Using MiniBatch and/or Sampling, Compute the final Q" << std::endl;
        }

        int molidx = 0, numvalmols = 0, numnonvalmols = 0;
        for (auto molData = molDataSet.begin(); molData != molDataSet.end(); ++molData, molidx++) {
            if (molData->GetGroup() == validation_group) {
                valQ += ComputeQ(*molData);
                numvalmols++;
            } else if (cfg->ga_minibatch_nth_size > 1 ||
                       sampling_method != USE_NO_SAMPLING) {
                Q += ComputeQ(*molData);
                numnonvalmols++;
            } else
                numnonvalmols++;
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master
        after = time(nullptr);
        std::string q_time_msg =
                "[M-Step] Finished Q compute: Time Elapsed = " +
                boost::lexical_cast<std::string>(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(q_time_msg.c_str());

        valQ = comm->collectQInMaster(valQ);
        if (cfg->ga_minibatch_nth_size > 1 ||
            sampling_method != USE_NO_SAMPLING) {
            if (comm->isMaster()) {
                Q += addRegularizersAndUpdateGradient(nullptr);
            }
            Q = comm->collectQInMaster(Q);
            Q = comm->broadcastQ(Q);
        }
        numvalmols = comm->collectSumInMaster(numvalmols);
        numnonvalmols = comm->collectSumInMaster(numnonvalmols);

        // Check for convergence
        double Qratio = fabs((Q - prevQ) / Q);
        double bestQRatio = fabs((Q - bestQ) / Q);
        if (comm->isMaster()) {
            std::string qdif_str = "[M-Step] Q_ratio= " + boost::lexical_cast<std::string>(Qratio) + " prev_Q=" +
                                   boost::lexical_cast<std::string>(prevQ) + "\n";
            qdif_str = "Best_Q_ratio= " + boost::lexical_cast<std::string>(bestQRatio) + " best_Q=" +
                       boost::lexical_cast<std::string>(bestQ) + "\n";

            qdif_str += "Q=" + boost::lexical_cast<std::string>(Q) +
                        " ValidationQ=" + boost::lexical_cast<std::string>(valQ) + " ";
            qdif_str +=
                    "Q_avg=" + boost::lexical_cast<std::string>(Q / numnonvalmols) +
                    " ValidationQ_avg=" + boost::lexical_cast<std::string>(valQ / numvalmols) + " ";
            writeStatus(qdif_str.c_str());
            comm->printToMasterOnly(qdif_str.c_str());
        }

        // check if make any progress yet
        // two conditions: 1. Qratio is less than 1e-15
        //                 2, Q has not improved compare to the best value so far
        const double ratio_cutoff = 1e-15;

        if (Qratio < ratio_cutoff || prevQ >= Q) {
            if (learning_rate > cfg->starting_step_size * 0.02) {
                learning_rate *= 0.5;
                count_no_progress = 0;
            } else {
                count_no_progress += 1;
            }
        } else {
            count_no_progress = 0;
        }

        // only save trhe best Q so far
        if (bestQ < Q) {
            bestQ = Q;
            // Write the params
            if (comm->isMaster()) {
                std::string progress_str = "[M-Step] Found Better Q: "
                                           + boost::lexical_cast<std::string>(bestQ) + " Write to File";
                comm->printToMasterOnly(progress_str.c_str());
                writeParamsToFile(iter_out_param_filename);
                writeParamsToFile(out_param_filename);
            }
        }

        prevQ = Q;
        // check if EM meet halt flag
        if (bestQRatio < cfg->em_converge_thresh || count_no_progress >= 3) {
            if (sampling_method != USE_NO_SAMPLING && cfg->reset_sampling) {
                if (comm->isMaster())
                    std::cout << "[Reset] Turn off sampling" << std::endl;
                sampling_method = USE_NO_SAMPLING;
                learning_rate = cfg->starting_step_size * cfg->reset_sampling_lr_ratio;
                count_no_progress = 0;
            } else {
                comm->printToMasterOnly(("EM Converged after " +
                                         boost::lexical_cast<std::string>(iter) +
                                         " iterations")
                                                .c_str());
                break;
            }
        }
        iter++;
    }


    return bestQ;
}