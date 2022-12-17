/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Config.cpp
#
# Description: 	Structs, functions, defaults for setting general 
#			    configuration data.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Config.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <boost/filesystem/operations.hpp>


void initDefaultConfig(config_t &cfg, char *argv_zero) {

    cfg.lambda = DEFAULT_LAMBDA;
    cfg.ga_method = DEFAULT_USE_ADAM_FOR_GA;
    cfg.em_converge_thresh = DEFAULT_EM_CONVERGE_THRESH;
    cfg.ga_converge_thresh = DEFAULT_GA_CONVERGE_THRESH;
    cfg.model_depth = DEFAULT_MODEL_DEPTH;
    cfg.abs_mass_tol = DEFAULT_ABS_MASS_TOL;
    cfg.ppm_mass_tol = DEFAULT_PPM_MASS_TOL;
    cfg.num_em_restarts = DEFAULT_NUM_EM_RESTARTS;
    cfg.starting_step_size = DEFAULT_LEARNING_RATE;
    cfg.ending_step_size = -DEFAULT_LEARNING_RATE;
    cfg.decay_rate = DEFAULT_DECAY_RATE;
    cfg.spectrum_depths.clear();
    cfg.spectrum_weights.clear();
    cfg.dv_spectrum_indexes.clear();
    cfg.fg_depth = DEFAULT_FRAGGRAPH_DEPTH;
    cfg.ionization_mode = DEFAULT_IONIZATION_MODE;
    cfg.update_bias_first = 0;
    cfg.param_init_type = PARAM_DEFAULT_INIT;
    cfg.use_lower_energy_params_for_init = 0;
    cfg.include_isotopes = DEFAULT_INCLUDE_ISOTOPES;
    cfg.isotope_thresh = DEFAULT_ISOTOPE_THRESH;
    cfg.allow_frag_detours = DEFAULT_ALLOW_FRAG_DETOURS;
    cfg.max_ring_breaks = DEFAULT_MAX_RING_BREAKS;
    cfg.theta_function = DEFAULT_THETA_FUNCTION;
    cfg.ga_minibatch_nth_size = DEFAULT_GA_MINIBATCH_NTH_SIZE;
    cfg.ga_max_iterations = DEFAULT_GA_MAX_ITERATIONS;
    cfg.ga_momentum = DEFAULT_GA_MOMENTUM;
    cfg.ga_adam_beta_1 = DEFAULT_ADAM_BETA_1;
    cfg.ga_adam_beta_2 = DEFAULT_ADAM_BETA_2;
    cfg.ga_adam_eps = DEFAULT_EPS;
    cfg.ga_adamw_w = DEFAULT_LAMBDA;
    cfg.ga_adadelta_rho = DEFAULT_ADADELTA_RHO;
    cfg.ga_decay_method = USE_NO_DECAY;
    cfg.exp_decay_k = DEFAULT_EXP_DECAY_K;
    cfg.step_decay_drop = DEFAULT_STEP_DECAY_DROP;
    cfg.step_decay_epochs_drop = DEFAULT_STEP_DECAY_EPOCHS_DROP;
    cfg.obs_function = DEFAULT_OBS_FUNCTION;
    cfg.include_h_losses = DEFAULT_INCLUDE_H_LOSSES;
    cfg.include_precursor_h_losses_only = DEFAULT_INCLUDE_PRECURSOR_H_LOSSES_ONLY;
    cfg.fragraph_compute_timeout_in_secs = DEFAULT_FRAGGRAPH_COMPUTE_TIMEOUT_IN_SECS;
    cfg.ga_use_best_q = DEFAULT_USE_BEST_Q_IN_GA;
    cfg.ga_sampling_method = USE_NO_SAMPLING;
    cfg.ga_sampling_method2 = USE_NO_SAMPLING;
    cfg.ga_reset_sampling = false;
    cfg.ga_sampling_max_selection = 100;
    cfg.ga_diff_sampling_peak_num = 30;
    cfg.ga_diff_sampling_difference = 0.1;
    cfg.disable_cross_val_metrics = false;
    cfg.disable_training_metrics = false;
    cfg.disable_cpu_usage_metrics = true;
    cfg.em_no_progress_count = 2;
    cfg.ga_no_progress_count = 3;
    cfg.collected_all_used_idx = false;
    cfg.em_max_iterations = DEFAULT_EM_MAX_ITERATIONS;
    cfg.allow_intermediate_peak = false;
    cfg.allow_cyclization = false;
    cfg.use_log_scale_peak = false;
    cfg.use_iterative_fg_gen = false;
    cfg.default_predicted_peak_min = 1;
    cfg.default_predicted_peak_max = 30;
    cfg.default_predicted_min_intensity = 0.0;
    cfg.default_postprocessing_energy = 80.0;
    cfg.default_mz_decimal_place = 5;
    cfg.intensity_msg_weight = 0.01;

    if (argv_zero == nullptr) {
        cfg.isotope_pattern_file = "ISOTOPE.DAT";
    }else{
        auto exeuc_dir = boost::filesystem::system_complete(argv_zero).remove_filename();
        cfg.isotope_pattern_file = exeuc_dir.append("ISOTOPE.DAT").string();
    }
}

void initConfig(config_t &cfg, std::string &filename,char *argv_zero, bool report_all) {

    std::string line, name;
    double value;
    std::ifstream ifs(filename.c_str(), std::ifstream::in);

    initDefaultConfig(cfg, argv_zero);

    //Read the config file into the paramater update config structure
    if (!ifs) std::cout << "Could not open file " << filename << std::endl;
    while (ifs.good()) {

        getline(ifs, line);
        if (line.size() < 3) continue;    //in case of empty line
        if ('#' == line[0]) continue; // in case comments in config
        std::stringstream ss1(line);
        ss1 >> name >> value;

        if (name == "lambda") cfg.lambda = value;
        else if (name == "ionization_mode") cfg.ionization_mode = (int) value;
        //else if (name == "converge_count_thresh") cfg.converge_count_thresh = (int) value;
        else if (name == "em_converge_thresh") cfg.em_converge_thresh = value;
        else if (name == "ga_converge_thresh") cfg.ga_converge_thresh = value;
        else if (name == "update_bias_first") cfg.update_bias_first = (int) value;
        else if (name == "model_depth") cfg.model_depth = (unsigned int) value;
        else if (name == "spectrum_depth") cfg.spectrum_depths.push_back((unsigned int) value);
        else if (name == "spectrum_weight") cfg.spectrum_weights.push_back((double) value);
        else if (name == "abs_mass_tol") cfg.abs_mass_tol = (double) value;
        else if (name == "ppm_mass_tol") cfg.ppm_mass_tol = (double) value;
        else if (name == "num_em_restarts") cfg.num_em_restarts = (int) value;
        else if (name == "starting_step_size") cfg.starting_step_size = (float) value;
        else if (name == "ending_step_size") cfg.ending_step_size = (float) value;
        else if (name == "decay_rate") cfg.decay_rate = (double) value;
        else if (name == "fg_depth") cfg.fg_depth = (int) value;
        else if (name == "allow_frag_detours") cfg.allow_frag_detours = (int) value;
        else if (name == "max_ring_breaks") cfg.max_ring_breaks = (int) value;
        else if (name == "include_isotopes") cfg.include_isotopes = (int) value;
        else if (name == "isotope_thresh") cfg.isotope_thresh = (double) value;
        else if (name == "ga_method") cfg.ga_method = (int) value;
        else if (name == "param_init_type") cfg.param_init_type = (int) value;
        else if (name == "use_lower_energy_params_for_init") cfg.use_lower_energy_params_for_init = (int) value;
        else if (name == "theta_function") cfg.theta_function = (int) value;
        else if (name == "theta_nn_hlayer_num_nodes") cfg.theta_nn_hlayer_num_nodes.push_back((int) value);
        else if (name == "theta_nn_layer_act_func_ids") cfg.theta_nn_layer_act_func_ids.push_back((int) value);
        else if (name == "nn_layer_dropout_probs") cfg.nn_layer_dropout_probs.push_back((float) value);
        else if (name == "nn_layer_freeze") cfg.nn_layer_is_frozen_flags.push_back((bool) value);
        else if (name == "ga_minibatch_nth_size") cfg.ga_minibatch_nth_size = (int) value;
        else if (name == "ga_max_iterations") cfg.ga_max_iterations = (int) value;
        else if (name == "ga_momentum") cfg.ga_momentum = (double) value;
        else if (name == "ga_adam_beta_1") cfg.ga_adam_beta_1 = (double) value;
        else if (name == "ga_adam_beta_2") cfg.ga_adam_beta_2 = (double) value;
        else if (name == "ga_adam_eps") cfg.ga_adam_eps = (double) value;
        else if (name == "ga_adamw_w") cfg.ga_adamw_w = (double) value;
        else if (name == "ga_adam_ga_adam_use_amsgrad") cfg.ga_adam_use_amsgrad = (bool) value;
        else if (name == "ga_adadelta_rho") cfg.ga_adadelta_rho = (double) value;
        else if (name == "ga_decay_method") cfg.ga_decay_method = (int) value;
        else if (name == "exp_decay_k") cfg.exp_decay_k = (double) value;
        else if (name == "step_decay_drop") cfg.step_decay_drop = (double) value;
        else if (name == "step_decay_epochs_drop") cfg.step_decay_epochs_drop = (double) value;
        else if (name == "obs_function") cfg.obs_function = (int) value;
        else if (name == "include_h_losses") cfg.include_h_losses = (int) value;
        else if (name == "include_precursor_h_losses_only") cfg.include_precursor_h_losses_only = (int) value;
        else if (name == "fragraph_compute_timeout_in_secs") cfg.fragraph_compute_timeout_in_secs = (int) value;
        else if (name == "ga_use_best_q") cfg.ga_use_best_q = (int) value;
        else if (name == "ga_sampling_method") cfg.ga_sampling_method = (int) value;
        else if (name == "ga_sampling_method2") cfg.ga_sampling_method2 = (int) value;
        else if (name == "ga_reset_sampling") cfg.ga_reset_sampling = (bool) value;
        else if (name == "ga_sampling_max_selection") cfg.ga_sampling_max_selection = (int) value;
        else if (name == "ga_diff_sampling_peak_num") cfg.ga_diff_sampling_peak_num = (int) value;
        else if (name == "ga_diff_sampling_difference") cfg.ga_diff_sampling_difference = (double) value;
        else if (name == "disable_cross_val_metrics") cfg.disable_cross_val_metrics = (int) value;
        else if (name == "disable_training_metrics") cfg.disable_training_metrics = (int) value;
        else if (name == "disable_cpu_usage_metrics") cfg.disable_cpu_usage_metrics = (bool) value;
        else if (name == "em_no_progress_count") cfg.em_no_progress_count = (int) value;
        else if (name == "ga_no_progress_count") cfg.ga_no_progress_count = (int) value;
        else if (name == "collected_all_used_idx") cfg.collected_all_used_idx = (bool) value;
        else if (name == "em_max_iterations") cfg.em_max_iterations = (int) value;
        else if (name == "allow_intermediate_peak") cfg.allow_intermediate_peak = (bool) value;
        else if (name == "allow_cyclization") cfg.allow_cyclization = (bool) value;
        else if (name == "use_log_scale_peak") cfg.use_log_scale_peak = (bool) value;
        else if (name == "use_iterative_fg_gen") cfg.use_iterative_fg_gen = (bool) value;
        else if (name == "default_predicted_peak_min") cfg.default_predicted_peak_min = (int) value;
        else if (name == "default_predicted_peak_max") cfg.default_predicted_peak_max = (int) value;
        else if (name == "default_predicted_min_intensity") cfg.default_predicted_min_intensity = (double) value;
        else if (name == "default_postprocessing_energy") cfg.default_postprocessing_energy = (double) value;
        else if (name == "default_mz_decimal_place") cfg.default_mz_decimal_place = (int) value;
        else if (name == "intensity_msg_weight") cfg.intensity_msg_weight = (double) value;
        else std::cout << "Warning: Unknown parameter configuration identifier " << name << std::endl;
    }
    ifs.close();

    if (cfg.spectrum_depths.size() != cfg.spectrum_weights.size())
        std::cout << "Warning: Mismatch between size of spectrum depths and weights" << std::endl;

    if (cfg.theta_function == NEURAL_NET_THETA_FUNCTION) {
        if (cfg.theta_nn_layer_act_func_ids.size() < cfg.theta_nn_hlayer_num_nodes.size() + 1) {
            std::cout
                    << "Warning: Activations function types not specified for all neural net layers, using default activations for unspecified layers."
                    << std::endl;
            while (cfg.theta_nn_layer_act_func_ids.size() < cfg.theta_nn_hlayer_num_nodes.size() + 1)
                cfg.theta_nn_layer_act_func_ids.push_back(DEFAULT_NN_ACTIVATION_FUNCTION);
        } else if (cfg.theta_nn_layer_act_func_ids.size() != cfg.theta_nn_hlayer_num_nodes.size() + 1) {
            std::cout << "Warning: More activations function types than neural net layers, ignoring some activations."
                      << std::endl;
            cfg.theta_nn_layer_act_func_ids.resize(cfg.theta_nn_hlayer_num_nodes.size() + 1);
        }
        for (int i = 0; i < cfg.theta_nn_hlayer_num_nodes.size(); i++) {
            if (cfg.theta_nn_layer_act_func_ids[i] == RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION &&
                cfg.theta_nn_hlayer_num_nodes[i] % 2 > 0) {
                std::cout << "Warning: Invalid to have odd number of nodes for ReLU hidden layer. Adding extra node."
                          << std::endl;
                cfg.theta_nn_hlayer_num_nodes[i] += 1;
            }
        }
        cfg.theta_nn_hlayer_num_nodes.push_back(1);    //Last layer has one node
    }
    // set default  ending_step_size
    if ( cfg.ending_step_size < 0)
        cfg.ending_step_size = cfg.starting_step_size * 0.25f;

    initDerivedConfig(cfg);

    //Report config parameters
    if (report_all) {
        if (cfg.ionization_mode == POSITIVE_ESI_IONIZATION_MODE)
            std::cout << "Positive ESI Ionization Mode" << std::endl;
        else if (cfg.ionization_mode == NEGATIVE_ESI_IONIZATION_MODE)
            std::cout << "Negative ESI Ionization Mode" << std::endl;
        else if (cfg.ionization_mode == POSITIVE_EI_IONIZATION_MODE)
            std::cout << "Positive EI Ionization Mode" << std::endl;
        else {
            std::cout << "Warning: Unknown Ionization Mode, reverting to default mode (positive)!" << std::endl;
            cfg.ionization_mode = DEFAULT_IONIZATION_MODE;
        }
        if (cfg.include_isotopes) {
            std::cout << "ISOTOPE.DAT file " << cfg.isotope_pattern_file << std::endl;
            std::cout << "Including fragment isotopes above intensity " << cfg.isotope_thresh << std::endl;
        }

        else std::cout << "Not including fragment isotopes" << std::endl;
        if (cfg.include_h_losses || cfg.include_precursor_h_losses_only) {
            std::cout << "Including Hydrogen losses";
            if (cfg.include_precursor_h_losses_only) std::cout << " from precursor only";
            std::cout << std::endl;
        }
        if (cfg.allow_intermediate_peak)
            std::cout << "Allowing intermediate peaks" << std::endl;
        if (cfg.allow_cyclization)
            std::cout << "Allowing cyclization" << std::endl;
        if (cfg.use_log_scale_peak)
            std::cout << "Using log scale peak" << std::endl;
        if (cfg.use_iterative_fg_gen)
            std::cout << "Using iterative fragmentation graph generation" << std::endl;

        std::cout << "Predicted peak num limited to [" << cfg.default_predicted_peak_min << ","
            << cfg.default_predicted_peak_max << "]" << std::endl;

        std::cout << "Predicted peak min intensity " << cfg.default_predicted_min_intensity << "%" << std::endl;
        std::cout << "Postprocessing energy " << cfg.default_postprocessing_energy << "%" << std::endl;
        std::cout << "Predicted peak mz to " << cfg.default_mz_decimal_place << " decimal place " << std::endl;

        if (cfg.param_init_type == PARAM_RANDOM_INIT) std::cout << "Using Random Parameter Initialisation" << std::endl;
        else if (cfg.param_init_type == PARAM_FULL_ZERO_INIT) std::cout << "Using Full Zero Initialisation" << std::endl;
        else if (cfg.param_init_type == PARAM_ZERO_INIT) std::cout << "Using Zero Initialisation (non-zero Bias)" << std::endl;
        else if (cfg.param_init_type == PARAM_NORMAL_INIT) std::cout << "Using NORMAL Initialisation FOR NN" << std::endl;
        else
            std::cout << "Warning: Unknown parameter initialization, revering to default mode (full random)!"
                      << std::endl;

        std::cout << "EM Restart Times: " << cfg.num_em_restarts << std::endl;
        std::cout << "Using EM Convergence Threshold " << cfg.em_converge_thresh << std::endl;
        std::cout << "Using Lambda " << cfg.lambda << std::endl;
        std::cout << "Using Intensity Weight " << cfg.intensity_msg_weight << " During Inference" << std::endl;


        if (USE_MOMENTUM_FOR_GA == cfg.ga_method) {
            std::cout << "Using simple gradient ascent implementation" << std::endl;
            std::cout << "Using Starting Step Size " << cfg.starting_step_size << " momentum " << cfg.ga_momentum
                      << " decay rate " << cfg.decay_rate << " Ending Step Size " << cfg.ending_step_size
                      << std::endl;
        } else if (USE_ADADELTA_FOR_GA == cfg.ga_method) {
            std::cout << "Using AdaDelta implementation" << std::endl;
            std::cout << "Using Starting Step Size " << cfg.starting_step_size << " decay rate " << cfg.decay_rate
                      << " eps " << cfg.ga_adam_eps << " Ending Step Size " << cfg.ending_step_size
                      << std::endl;
        } else if (USE_ADAM_FOR_GA == cfg.ga_method
                   || USE_ADABELIEF_FOR_GA == cfg.ga_method
                   || USE_ADAMW_FOR_GA == cfg.ga_method) {
            if (USE_ADAM_FOR_GA == cfg.ga_method)
                std::cout << "Using Adam implementation" << std::endl;
            else if (USE_ADABELIEF_FOR_GA == cfg.ga_method)
                std::cout << "Using AdaBelief implementation" << std::endl;
            else if (USE_ADAMW_FOR_GA == cfg.ga_method){
                std::cout << "Using AdamW implementation" << std::endl;
                if(cfg.lambda != 0.0){
                    cfg.lambda = 0.0;
                    std::cout << "Set lambda to " << cfg.lambda << ", adamW use Weight Decay not L2" << std::endl;
                }
            }
            std::cout << "Using starting step size " << cfg.starting_step_size  << " ending step size " << cfg.ending_step_size
                     << " beta1  " << cfg.ga_adam_beta_1 << " beta2 " << cfg.ga_adam_beta_2 
                     << " eps " << cfg.ga_adam_eps << " amsgrad " << cfg.ga_adam_use_amsgrad;
            if (USE_ADAMW_FOR_GA == cfg.ga_method)
                std::cout << " ga_adamw_w " << cfg.ga_adamw_w;
            std::cout << std::endl;
        }
        std::cout << "Using GA max iterations " << cfg.ga_max_iterations << std::endl;
        std::cout << "Using GA Convergence Threshold " << cfg.ga_converge_thresh << std::endl;
        std::cout << "Using GA mini batch taking 1 in " << cfg.ga_minibatch_nth_size << " of processor data"
                  << std::endl;

        switch (cfg.ga_decay_method) {
            case USE_DEFAULT_DECAY:
                std::cout << "Using Time Based Learning Rate Decay Method. decay rate " << cfg.decay_rate << std::endl;
                break;
            case USE_EXP_DECAY:
                std::cout << "Using Exponential Learning Rate Decay Method. k: " << cfg.exp_decay_k << std::endl;
                break;
            case USE_STEP_DECAY:
                std::cout << "Using Step Decay Learning Rate Method. drop: " << cfg.step_decay_drop << " epochs drop: "
                          << cfg.step_decay_epochs_drop << std::endl;
                break;
            case USE_NO_DECAY:
            default:
                std::cout << "NOT Using Decay Method" << std::endl;
        }
        if (cfg.ga_use_best_q)
            std::cout << "Using Best Q instead of Prev Q in GA" << std::endl;
        std::cout << "Using Fragmentation Graph Depth " << cfg.fg_depth << std::endl;

        std::cout << "Using ";
        printSamplingConfig(cfg.ga_sampling_method, cfg);
        if (cfg.ga_reset_sampling) {
            std::cout << "Reset Sampling Enabled , Reset to ";
            printSamplingConfig(cfg.ga_sampling_method2, cfg);
        }

        if (cfg.allow_frag_detours) {
            std::cout << "Allowing fragmentation detours " << std::endl;
        } else {
            std::cout << "Disallowing fragmentation detours ";
        }
        std::cout << "Maximum Ring Breaks " << cfg.max_ring_breaks << std::endl;
        std::cout << "Using Model Depth " << cfg.model_depth << std::endl;
        std::cout << "Using Spectrum Depths and Weights: ";
        for (unsigned int i = 0; i < cfg.spectrum_depths.size(); i++)
            std::cout << "(" << cfg.spectrum_depths[i] << "," << cfg.spectrum_weights[i] << ") ";
        std::cout << std::endl;
        std::cout << "Using Absolute mass tolerance " << cfg.abs_mass_tol << std::endl;
        std::cout << "Using PPM mass tolerance " << cfg.ppm_mass_tol << std::endl;
        if (cfg.use_lower_energy_params_for_init)
            std::cout << "Initialising higher energy params with those of one level lower" << std::endl;
        if (cfg.theta_function == LINEAR_THETA_FUNCTION) std::cout << "Using linear function for theta" << std::endl;
        else if (cfg.theta_function == NEURAL_NET_THETA_FUNCTION) {
            std::cout << "Using neural net for theta with " << cfg.theta_nn_hlayer_num_nodes.size()
                      << " hidden layers: " << std::endl;

            std::cout << "Hidden layer nodes:   ";        
            for (int i = 0; i < cfg.theta_nn_hlayer_num_nodes.size(); i++)
                std::cout << std::setw(10) <<  cfg.theta_nn_hlayer_num_nodes[i];
            std::cout << std::endl;
            std::cout << "Activation functions: ";
            for (int i = 0; i < cfg.theta_nn_layer_act_func_ids.size(); i++){
                switch (cfg.theta_nn_layer_act_func_ids[i]){
                    case RELU_NN_ACTIVATION_FUNCTION:
                        std::cout << std::setw(10) << "Relu";
                        break;
                    case LEAKY_RELU_NN_ACTIVATION_FUNCTION:
                        std::cout << std::setw(10) << "Leaky Relu";
                        break;
                    case RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION: {
                        if (i % 2 == 0)
                            std::cout << std::setw(10) << "Neg Relu";
                        else
                            std::cout << std::setw(10) << "Relu ";
                        break;
                    }
                    default:
                        std::cout << std::setw(10) << "Linear";
                }
            }
            std::cout << std::endl;

            std::cout << "Hidden layer dropout: ";         
            for (int i = 0; i < cfg.nn_layer_dropout_probs.size(); i++)
                std::cout << std::setw(10) <<  cfg.nn_layer_dropout_probs[i];
            std::cout << std::endl;

            std::cout << "Hidden layer freeze:  ";         
            for (int i = 0; i < cfg.nn_layer_is_frozen_flags.size(); i++)
                std::cout << std::setw(10) <<  cfg.nn_layer_is_frozen_flags[i];
            std::cout << std::endl;
        }
        if (cfg.obs_function == NORMAL_OBS_FUNCTION)
            std::cout << "Using normally distributed observation function" << std::endl;
        else if (cfg.obs_function == UNIFORM_OBS_FUNCTION)
            std::cout << "Using windowed uniform distribution observation function" << std::endl;
        else {
            std::cout << "Warning: Unrecognised observation function (" << cfg.obs_function << "). Using default"
                      << std::endl;
            cfg.obs_function = DEFAULT_OBS_FUNCTION;
        }
        if (cfg.fragraph_compute_timeout_in_secs > 0)
            std::cout << "Timeout set on fragment graph computation to " << cfg.fragraph_compute_timeout_in_secs
                      << " seconds" << std::endl;
        if (cfg.use_log_scale_peak)
            std::cout << "Use log scale peaks" << std::endl;
    }
}

void printSamplingConfig(const int &sampling_method, config_t &cfg) {
    switch (sampling_method) {
        case USE_RANDOM_SAMPLING:
            std::cout << "Random Sampling on transitions  with max selection cut off at "
                      << cfg.ga_sampling_max_selection << std::endl;
            break;
        case USE_GRAPH_RANDOM_WALK_SAMPLING:
            std::cout << "Graph Random Walk Sampling with "
                      << cfg.ga_sampling_max_selection << " iterations" << std::endl;
            break;
        case USE_DIFFERENCE_SAMPLING_BFS_CO:
            std::cout << "Difference sampling BFS Child Only with max selected peak count at "
                      << cfg.ga_diff_sampling_peak_num << std::endl;
            break;
        case USE_DIFFERENCE_SAMPLING_BFS:
            std::cout << "Difference sampling BFS with max selected peak count at "
                      << cfg.ga_diff_sampling_peak_num << std::endl;
            break;
        case USE_NO_SAMPLING:
        default:
            std::cout << "NO Sampling Method" << std::endl;
    }
}

void initDerivedConfig(config_t &cfg, int se_energy) {

    //Derived Parameters
    cfg.map_d_to_energy.resize(cfg.model_depth);
    int energy = 0;
    if (se_energy > 0) energy = se_energy;
    for (unsigned int d = 0; d < cfg.model_depth; d++) {
        std::vector<int>::iterator it = cfg.spectrum_depths.begin();
        for (; it != cfg.spectrum_depths.end(); ++it) {
            if (d == *it) energy++;
        }
        cfg.map_d_to_energy[d] = energy;
    }

    cfg.dv_spectrum_depths = cfg.spectrum_depths;
    cfg.dv_spectrum_weights = cfg.spectrum_weights;

    //Re-normalise weights
    double sum = 0.0;
    std::vector<double>::iterator it = cfg.dv_spectrum_weights.begin();
    for (; it != cfg.dv_spectrum_weights.end(); ++it) sum += *it;
    it = cfg.dv_spectrum_weights.begin();
    for (; it != cfg.dv_spectrum_weights.end(); ++it) *it = *it / sum;

    //Set spectrum indexes
    cfg.dv_spectrum_indexes.clear();
    if (se_energy >= 0) cfg.dv_spectrum_indexes.push_back(se_energy);
    else {
        for (int i = 0; i < cfg.dv_spectrum_depths.size(); i++)
            cfg.dv_spectrum_indexes.push_back(i);
    }

}

void initSingleEnergyConfig(config_t &se_cfg, config_t &cfg, int energy) {

    se_cfg = cfg;

    //Adjust the depth parameters to include only one spectrum
    se_cfg.model_depth = cfg.spectrum_depths[energy];
    se_cfg.spectrum_depths.resize(1);
    se_cfg.spectrum_depths[0] = cfg.spectrum_depths[energy];
    se_cfg.spectrum_weights.resize(1);
    se_cfg.spectrum_weights[0] = 1.0;

    //Re-derive the derived parameters
    initDerivedConfig(se_cfg, energy);
}