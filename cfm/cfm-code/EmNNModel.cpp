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

#include "mpi.h"
#include "EmNNModel.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>


EmNNModel::EmNNModel(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename,
             std::string initial_params_filename) : EmModel(a_cfg, an_fc, a_status_filename, initial_params_filename) {
    int num_energies_to_include = cfg->map_d_to_energy.back() + 1;
    if (!initial_params_provided) {
        nn_param = boost::shared_ptr<NNParam>(
                new NNParam(fc->getFeatureNames(), num_energies_to_include, a_cfg->theta_nn_hlayer_num_nodes,
                            a_cfg->theta_nn_layer_act_func_ids));
        comm->printToMasterOnly("EM_NN: No initial params provided");
    } else {
        nn_param = boost::shared_ptr<NNParam>(new NNParam(initial_params_filename));
        while (nn_param->getNumEnergyLevels() < num_energies_to_include)
            nn_param->appendRepeatedPrevEnergyParams();
        std::string msg = "EM_NN: Initial params provided from " + initial_params_filename;
        comm->printToMasterOnly(msg.c_str());
    }
    param = nn_param;
}

//These functions are the same as for EM, but use NNParam, so need to be redefined here
void EmNNModel::initParams() {
    switch (cfg->param_init_type){
        case PARAM_FULL_ZERO_INIT:
            nn_param->fullZeroInit();
            break;
        case PARAM_ZERO_INIT:
            nn_param->zeroInit();
            break;
        case PARAM_NORMAL_INIT:
            nn_param->randomNormalInit();
            break;
        case NN_PARAM_VAR_SCALING_INIT:
            nn_param->varianceScalingInitializer();
            break;
        case PARAM_RANDOM_INIT:
        default:
            nn_param->randomUniformInit();
    }
}

void EmNNModel::computeThetas(MolData *moldata) {
    moldata->computeTransitionThetas(*nn_param);
}

void EmNNModel::writeParamsToFile(std::string &filename) {
    nn_param->saveToFile(filename);
}

//Gradient Computation using Backpropagation
double EmNNModel::computeAndAccumulateGradient(double *grads, int molidx, MolData &moldata, suft_counts_t &suft,
                                               bool record_used_idxs_only, std::set<unsigned int> &used_idxs,
                                               int sampling_method, double sampling_explore_rate) {

    double Q = 0.0;
    //const FragmentGraph *fg = moldata.getFragmentGraph();
    unsigned int num_transitions = moldata.getNumTransitions();
    unsigned int num_fragments = moldata.getNumFragments();

    int offset = num_transitions;
    unsigned int layer2_offset = nn_param->getSecondLayerWeightOffset();
    int weights_per_energy = nn_param->getNumWeightsPerEnergyLevel();

    if (!moldata.hasComputedGraph()) return Q;

    suft_t *suft_values = &(suft.values[molidx]);

    //Collect energies to compute
    std::vector<unsigned int> energies;
    unsigned int energy;
    int prev_energy = -1;
    for (unsigned int d = 0; d < cfg->model_depth; d++) {
        energy = cfg->map_d_to_energy[d];
        if (energy != prev_energy) energies.push_back(energy);
        prev_energy = energy;
    }

    //Compute the gradients
    std::vector<unsigned int>::iterator eit = energies.begin();
    for (; eit != energies.end(); ++eit) {
        energy = *eit;

        unsigned int grad_offset = energy * nn_param->getNumWeightsPerEnergyLevel();
        unsigned int suft_offset = energy * (num_transitions + num_fragments);


        std::set<int> selected_trans_id;
        if(!record_used_idxs_only && sampling_method != USE_NO_SAMPLING)
            getRandomWalkedTransitions(moldata, sampling_method, sampling_explore_rate, energy, selected_trans_id);

        //Iterate over from_id (i)
        const tmap_t *from_map = moldata.getFromIdTMap();
        for (int from_idx = 0; from_idx < num_fragments; from_idx++) {

            const std::vector<int> *from_id_map = &(*from_map)[from_idx];
            std::vector<int> sampled_trans_id;

            if (sampling_method != USE_NO_SAMPLING && !record_used_idxs_only) {
                for (const auto &trans_id : (*from_map)[from_idx]) {
                    if (selected_trans_id.find(trans_id) != selected_trans_id.end()) {
                        sampled_trans_id.push_back(trans_id);
                    }
                }
                from_id_map = &sampled_trans_id;
            }

            if (from_id_map->empty())
                continue;

            unsigned int num_trans_from_id = from_id_map->size();
            std::vector<azd_vals_t> a_values(num_trans_from_id);
            std::vector<azd_vals_t> z_values(num_trans_from_id);

            //Compute the forward values, storing intermediate a and z values, and the combined denom of the rho term, and Q
            double denom = 1.0, nu_sum = 0.0;
            std::vector<int>::const_iterator it = from_id_map->begin();
            std::vector<const FeatureVector *> fvs(num_trans_from_id);
            std::vector<double> nu_terms(num_trans_from_id + 1);
            for (int idx = 0; it != from_id_map->end(); ++it, idx++) {
                fvs[idx] = moldata.getFeatureVectorForIdx(*it);
                double theta = nn_param->computeTheta(*fvs[idx], energy, z_values[idx], a_values[idx]);
                denom += exp(theta);
                nu_terms[idx] = (*suft_values)[*it + suft_offset];
                nu_sum += nu_terms[idx];
                Q += nu_terms[idx] * theta;
            }
            nu_terms[num_trans_from_id] = (*suft_values)[offset + from_idx + suft_offset];
            nu_sum += nu_terms[num_trans_from_id];
            Q -= log(denom) * nu_sum;

            //Compute the deltas
            std::vector<azd_vals_t> deltasA, deltasB;
            nn_param->computeDeltas(deltasA, deltasB, z_values, a_values, denom, energy);

            //Compute the unweighted gradients
            std::vector<std::vector<double> > unweighted_grads;
            std::set<unsigned int> from_id_used_idxs;
            nn_param->computeUnweightedGradients(unweighted_grads, from_id_used_idxs, fvs, deltasA, deltasB, a_values);

            //Accumulate the weighted gradients
            std::set<unsigned int>::iterator sit;
            it = from_id_map->begin();
            for (int idx = 0; idx <= num_trans_from_id; idx++) {
                //First layer
                double nu = nu_terms[idx];
                for (sit = from_id_used_idxs.begin(); sit != from_id_used_idxs.end(); ++sit)
                    *(grads + grad_offset + *sit) += nu * unweighted_grads[idx][*sit];
                //Remaining layers
                for (int i = layer2_offset; i < weights_per_energy; i++)
                    *(grads + grad_offset + i) += nu * unweighted_grads[idx][i];
            }

            if (record_used_idxs_only) {
                //Combine the used indexes for the first layer
                for (sit = from_id_used_idxs.begin(); sit != from_id_used_idxs.end(); ++sit)
                    used_idxs.insert(grad_offset + *sit);
            }
        }

        if (record_used_idxs_only) {
            //Add used indexes for the subsequent layers (all of them)
            for (int i = layer2_offset; i < weights_per_energy; i++)
                used_idxs.insert(grad_offset + i);
        }
    }
    return Q;
}

double EmNNModel::computeQ(int molidx, MolData &moldata, suft_counts_t &suft) {

    double Q = 0.0;
    //const FragmentGraph *fg = moldata.getFragmentGraph();
    unsigned int num_transitions = moldata.getNumTransitions();
    unsigned int num_fragments = moldata.getNumFragments();

    int offset = num_transitions;
    unsigned int layer2_offset = nn_param->getSecondLayerWeightOffset();
    int weights_per_energy = nn_param->getNumWeightsPerEnergyLevel();

    if (!moldata.hasComputedGraph()) return Q;

    suft_t *suft_values = &(suft.values[molidx]);

    //Collect energies to compute
    std::vector<unsigned int> energies;
    unsigned int energy;
    int prev_energy = -1;
    for (unsigned int d = 0; d < cfg->model_depth; d++) {
        energy = cfg->map_d_to_energy[d];
        if (energy != prev_energy) energies.push_back(energy);
        prev_energy = energy;
    }

    //Compute Q
    std::vector<unsigned int>::iterator eit = energies.begin();
    for (; eit != energies.end(); ++eit) {
        energy = *eit;

        unsigned int suft_offset = energy * (num_transitions + num_fragments);

        //Iterate over from_id (i)
        const tmap_t *from_map = moldata.getFromIdTMap();
        for (int from_idx = 0; from_idx < num_fragments; from_idx++) {

            const std::vector<int> *from_id_map = &(*from_map)[from_idx];
            unsigned int num_trans_from_id = from_id_map->size();

            //Compute the forward values, and the combined denom of the rho term, and Q
            double denom = 1.0, nu_sum = 0.0;
            std::vector<int>::const_iterator it = from_id_map->begin();
            std::vector<const FeatureVector *> fvs(num_trans_from_id);
            std::vector<double> nu_terms(num_trans_from_id + 1);
            for (int idx = 0; it != from_id_map->end(); ++it, idx++) {
                fvs[idx] = moldata.getFeatureVectorForIdx(*it);
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
    }
    return Q;
}


double EmNNModel::addRegularizersAndUpdateGradient(double *grads) {

    double Q = 0.0;
    auto it = ((MasterComms *) comm)->master_used_idxs.begin();
    for (; it != ((MasterComms *) comm)->master_used_idxs.end(); ++it) {
        double weight = nn_param->getWeightAtIdx(*it);
        Q -= 0.5 * cfg->lambda * weight * weight;
        if (grads != nullptr)
            *(grads + *it) -= cfg->lambda * weight;
    }

    //Remove part of the Bias terms (regularize the bias terms 100 times less than the other weights!)
    unsigned int weights_per_energy = nn_param->getNumWeightsPerEnergyLevel();
    std::vector<unsigned int> bias_indexes;
    nn_param->getBiasIndexes(bias_indexes);
    for (unsigned int energy = 0; energy < nn_param->getNumEnergyLevels(); energy++) {
        auto it = bias_indexes.begin();
        int offset = energy * weights_per_energy;
        for (; it != bias_indexes.end(); ++it) {
            double bias = nn_param->getWeightAtIdx(offset + *it);
            Q += 0.99 * 0.5 * cfg->lambda * bias * bias;
            if (grads != nullptr)
                *(grads + offset + *it) += 0.99 * cfg->lambda * bias;
        }
    }
    return Q;
}

