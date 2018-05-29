//
// Created by feiw on 26/05/18.
//

#include "DirectModel.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/exponential.hpp>

void DirectModel::computeAndAccumulateGradient(double *grads, MolData &moldata,
                                               bool record_used_idxs_only, std::set<unsigned int> &used_idxs,
                                               int sampling_method) {
    if (!moldata.hasComputedGraph())
        return;

    moldata.computeLogTransitionProbabilities();
	const FragmentGraph *fg = moldata.getFragmentGraph();
	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();

    // Collect energies to compute
    std::vector<unsigned int> energies;
    getEnergiesLevels(energies);

    // Compute the gradients
    for (auto energy_level : energies) {
		unsigned int grad_offset = energy_level * param->getNumWeightsPerEnergyLevel();

        for(const auto & peak : *moldata.getSpectrum(energy_level)){
            // get intensity or height of the spectrum
            double intensity = peak.intensity;
            // get possible path to this peak
            // anything out of mass tol range will not be considered
            std::vector<Path> paths;
            double mass_tol = getMassTol(cfg->abs_mass_tol,cfg->ppm_mass_tol,peak.mass);
            moldata.getPathes(paths, peak.mass, mass_tol);

            int num_transitions = moldata.getFragmentGraph()->getNumTransitions();
            for(auto & path: paths){
                // mean is the peak mass
                // std is the mass tol
                boost::math::normal_distribution<double> normal_dist(peak.mass, mass_tol);
                for(auto & trans_id : *path.getTransIds()){
                    if(trans_id < num_transitions){
						const FeatureVector *fv = moldata.getFeatureVectorForIdx(trans_id);
						for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
							auto fv_idx = *fv_it;

						}
                    } else {
                        //q += moldata.getLogPersistenceProbForIdx(energy_level, trans_id);
						const FeatureVector *fv = moldata.getFeatureVectorForIdx(trans_id);
						//double thete = moldata.getThetaForIdx(trans_id);
						for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
							auto fv_idx = *fv_it;

						}
                    }
                }
                // add observation term
                // log(boost::math::pdf(normal_dist, path.getDstMass()));
            }

        }
    }

    // Compute the latest transition thetas
    if (!record_used_idxs_only)
        moldata.computeNormalizedTransitionThetas(*param);
}

double DirectModel::computeQ(MolData &moldata) {

    double q = 0.0;
    if (!moldata.hasComputedGraph())
        return q;

    moldata.computeLogTransitionProbabilities();

    // Collect energies to compute
    std::vector<unsigned int> energies;
    getEnergiesLevels(energies);

    // Compute the gradients
    for (auto energy_level : energies) {
        for(const auto & peak : *moldata.getSpectrum(energy_level)) {
            // get intensity or height of the spectrum
            double intensity = peak.intensity;
            // get possible path to this peak
            // anything out of mass tol range will not be considered
            std::vector<Path> paths;
            double mass_tol = getMassTol(cfg->abs_mass_tol, cfg->ppm_mass_tol, peak.mass);
            moldata.getPathes(paths, peak.mass, mass_tol);

            int num_transitions = moldata.getFragmentGraph()->getNumTransitions();
            for (auto &path: paths) {
                // mean is the peak mass
                // std is the mass tol
                boost::math::normal_distribution<double> normal_dist(peak.mass, mass_tol);
                for (auto &trans_id : *path.getTransIds()) {
                    if (trans_id < num_transitions) {
                        q += moldata.getLogTransitionProbForIdx(energy_level, trans_id);
                    } else {
                        q += moldata.getLogPersistenceProbForIdx(energy_level, trans_id);
                    }
                }
                // add observation term
                q = log(boost::math::pdf(normal_dist, path.getDstMass()));
            }
            q *= intensity;
        }
    }
    // return negative log likehood
    return -q;
}

double DirectModel::addRegularizersAndUpdateGradient(double *grads) {
    return 0.0;
}

double DirectModel::updateParametersGradientAscent(std::vector<MolData> &data,
                                                   double learning_rate, int sampling_method) {
    return 0.0;
}

double DirectModel::trainModel(std::vector<MolData> &molDataSet, int group, std::string &out_param_filename) {
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
            if (mol.getGroup() != validation_group)
                mol.addNoise(cfg->noise_max, cfg->noise_sum, cfg->abs_mass_tol, cfg->ppm_mass_tol);
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
            if (molData->getGroup() == validation_group) {
                valQ += computeQ(*molData);
                numvalmols++;
            } else if (cfg->ga_minibatch_nth_size > 1 ||
                       sampling_method != USE_NO_SAMPLING) {
                Q += computeQ(*molData);
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