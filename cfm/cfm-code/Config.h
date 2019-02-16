/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Config.h
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

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <vector>
#include <string>

//Default Values

//Fragmentation depth in the initial fragmentation graph
static const int DEFAULT_FRAGGRAPH_DEPTH = 2;

//Fragmentation depth in the model
static const int DEFAULT_MODEL_DEPTH = 6;
static const int DEFAULT_MED_DEPTH = 4;
static const int DEFAULT_LOW_DEPTH = 2;

//Observation variance (parameters to the peak gaussian distribution -
// indicative of mass spec mass accuracy)
static const double DEFAULT_ABS_MASS_TOL = 0.01;
static const double DEFAULT_PPM_MASS_TOL = 10.0;

//Threshold below which the gradient needs to drop to end gradient ascent (in UpdateParameters)
static const double DEFAULT_GA_CONVERGE_THRESH = 0.0001;

//Regularization Constant
static const double DEFAULT_LAMBDA = 0.01;

// GA methods
static const int USE_SGD_FOR_GA = 0;
static const int USE_MOMENTUM_FOR_GA = 1;
static const int USE_ADAM_FOR_GA = 2;
static const int USE_ADAMW_FOR_GA = 3;
static const int USE_ADADELTA_FOR_GA = 4;
static const int USE_AMSGRAD_FOR_GA = 5;

//Use ADAM to do gradient ascent
static const int DEFAULT_USE_ADAM_FOR_GA = USE_ADAM_FOR_GA;

// LR Decay methods
static const int USE_NO_DECAY = 0;
static const int USE_DEFAULT_DECAY = 1;
static const int USE_EXP_DECAY = 2;
static const int USE_STEP_DECAY = 3;

// Decay
static const double DEFAULT_DECAY_RATE = 0.0;
static const double DEFAULT_EXP_DECAY_K = 0.1;
static const double DEFAULT_STEP_DECAY_DROP = 0.5;
static const double DEFAULT_STEP_DECAY_EPOCHS_DROP = 10;

// Learning Rate
static const double DEFAULT_LEARNING_RATE = 0.001;

//Momentum term used in Gradient Ascent
static const double DEFAULT_MOMENTUM_ALPHA = 0.1;

//ADAM
static const double DEFAULT_ADAM_BETA_1 = 0.9;
static const double DEFAULT_ADAM_BETA_2 = 0.999;

//AdaDelta
static const double DEFAULT_ADADELTA_LEARNING_RATE = 1.0;
static const double DEFAULT_ADADELTA_RHO = 0.95;

// EPS
static const double DEFAULT_EPS = 1e-8;

//Whether or not gradient ascent should be stochastic
static const int DEFAULT_STOCHASTIC_GA = 0;

//Max by which Q can change between iterations to call convergence 
static const double DEFAULT_EM_CONVERGE_THRESH = 0.001;

//Number of iterations gradient ascent must remain converged
//before being considered converged
static const int DEFAULT_CONVERGE_COUNT_THRESH = 1;

//Which IPFP-like algorithm to run during E-step
static const int DEFAULT_IPFP_ALGORITHM = 2;    //IPFP_WITH_OSC_ADJUST

//Threshold for testing convergence of the IPFP algorithm
static const double DEFAULT_IPFP_CONVERGE_THRESH = 0.005;

//Threshold for testing oscillatory convergence of the IPFP algorithm
static const double DEFAULT_IPFP_OSC_CONVERGE_THRESH = 0.999;

//Number of times to do a random restart of EM
static const int DEFAULT_NUM_EM_RESTARTS = 3;

//Default weight for undetected mass in the spectrum
static const double DEFAULT_LOST_PROB = 0.0;

//Default weight given to the next highest spectrum reading for unobserved fragments
static const double DEFAULT_BETWEEN_LOST_PROB = 1.0;

//Gradient Ascent Line Search Parameters
static const double DEFAULT_LINE_SEARCH_ALPHA = 0.1;
static const double DEFAULT_LINE_SEARCH_BETA = 0.5;
static const int DEFAULT_MAX_SEARCH_COUNT = 20;

static const int POSITIVE_ESI_IONIZATION_MODE = 1;
static const int NEGATIVE_ESI_IONIZATION_MODE = 2;
static const int POSITIVE_EI_IONIZATION_MODE = 3;

static const int DEFAULT_IONIZATION_MODE = POSITIVE_ESI_IONIZATION_MODE;
static const int DEFAULT_INCLUDE_ISOTOPES = 0;
static const double DEFAULT_ISOTOPE_THRESH = 1.0;    //On a normalized scale where the highest peak has height 100.0
static const int DEFAULT_INCLUDE_H_LOSSES = 1;        //For backwards compatibility with old ESI-MS/MS models
static const int DEFAULT_INCLUDE_PRECURSOR_H_LOSSES_ONLY = 0;

static const int PARAM_FULL_ZERO_INIT = 1;
static const int PARAM_ZERO_INIT = 2;
static const int PARAM_RANDOM_INIT = 3;
static const int PARAM_NORMAL_INIT = 4;
static const int NN_PARAM_VAR_SCALING_INIT = 5;

static const int PARAM_DEFAULT_INIT = PARAM_RANDOM_INIT;

static const int DEFAULT_ALLOW_FRAG_DETOURS = 1;
static const int DEFAULT_DO_PRELIM_BFS = 1;
static const int DEFAULT_MAX_RING_BREAKS = 2;

static const int MAX_BREAKABLE_RING_SIZE = 8;

//Mode for writing spectra to output
static const int NO_OUTPUT_MODE = 0;
static const int MSP_OUTPUT_MODE = 1;
static const int MGF_OUTPUT_MODE = 2;

//Maximum additional electron pairs to either side during fragmentation
static const int MAX_E_MOVE = 4;

static const int LINEAR_THETA_FUNCTION = 1;
static const int NEURAL_NET_THETA_FUNCTION = 2;
static const int DEFAULT_THETA_FUNCTION = LINEAR_THETA_FUNCTION;

static const int LINEAR_NN_ACTIVATION_FUNCTION = 0;
static const int RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION= 1;
static const int RELU_NN_ACTIVATION_FUNCTION = 2;
static const int LEAKY_RELU_NN_ACTIVATION_FUNCTION = 3;

static const int DEFAULT_NN_ACTIVATION_FUNCTION = RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION;

static const double DEFAULT_GA_MOMENTUM = 0.9;

static const int DEFAULT_GA_MINIBATCH_NTH_SIZE = 1;
static const int DEFAULT_GA_MAX_ITERATIONS = 20;

static const int NORMAL_OBS_FUNCTION = 1;
static const int UNIFORM_OBS_FUNCTION = 2;
static const int DEFAULT_OBS_FUNCTION = NORMAL_OBS_FUNCTION;

//Timeout settings
static const int DEFAULT_FRAGGRAPH_COMPUTE_TIMEOUT_IN_SECS = -1;    //-1 is no timeout.

//Random Sample settings
static const int DEFAULT_NOT_USE_GRAPH_PRUNING = 0;
static const int DEFAULT_USE_BEST_Q_IN_GA = 0;
static const int USE_NO_SAMPLING = 0;
static const int USE_GRAPH_RANDOM_WALK_SAMPLING = 1;
static const int USE_GRAPH_WEIGHTED_RANDOM_WALK_SAMPLING = 2;
static const int USE_DIFFERENCE_SAMPLING = 3;
static const double DEFAULT_GRAPH_SAMPLING_K = 1.0;

//Configuration
struct config_t {

    //Fragment Graph Configuration
    int fg_depth;
    int allow_frag_detours;
    int max_ring_breaks;
    int include_h_losses;
    int include_precursor_h_losses_only;

    //Use Single Energy CFM (rather than Combined Energy)
    int ionization_mode;
    int include_isotopes;
    double isotope_thresh;

    // For sampling
    int use_graph_pruning;
    bool reset_sampling;
    double reset_sampling_lr_ratio;
    int ga_sampling_method;
    double ga_graph_sampling_k;
    double ga_sampling_explore_weight;
    int ga_diff_sampling_peak_num;
    double ga_diff_sampling_difference;
    //int max_transitions_per_iteration;

    //Model Level Configuration
    unsigned int model_depth; //Total Depth
    std::vector<int> spectrum_depths;
    std::vector<double> spectrum_weights;
    double abs_mass_tol;
    double ppm_mass_tol;
    std::vector<int> map_d_to_energy;    //Derived parameter
    std::vector<int> dv_spectrum_depths;    //Either a direct copy, or interpolated values.
    std::vector<int> dv_spectrum_indexes;    //Index of each spectrum in the list of spectra in MolData
    std::vector<double> dv_spectrum_weights;
    int obs_function;    //Function used for the observation function

    //Theta function Configuration (Linear or Neural Net)
    int theta_function;
    //Neural Net Configuration
    std::vector<int> theta_nn_hlayer_num_nodes;
    std::vector<int> theta_nn_layer_act_func_ids;
    std::vector<double> nn_layer_dropout_probs;

    //EM Configuration
    double em_converge_thresh;
    int num_em_restarts;
    int use_lower_energy_params_for_init;
    int param_init_type;

    //Gradient Ascent Configuration
    double lambda;    //Regularization n
    int ga_method;
    int ga_decay_method;
    double ga_converge_thresh;

    double starting_step_size;
    int ga_max_iterations;
    int ga_use_best_q;

    // Decay
    double decay_rate;
    double exp_decay_k;
    double step_decay_drop;
    double step_decay_epochs_drop;
    // For Momentum
    double ga_momentum;
    // For ADAM
    double ga_adam_beta_1;
    double ga_adam_beta_2;
    double ga_adam_eps;
    // For ADAMW
    double ga_adamw_w;

    // For adadelta
    double ga_adadelta_rho;

    int max_search_count;
    int update_bias_first;
    int ga_minibatch_nth_size;
    double ga_dropout_delta;

    int fragraph_compute_timeout_in_secs;

    bool disable_cross_val_metrics;
    bool disable_training_metrics;
};

void initDefaultConfig(config_t &cfg);

void initConfig(config_t &cfg, std::string &filename, bool report_all = false);

void initDerivedConfig(config_t &cfg, int energy = -1);

void initSingleEnergyConfig(config_t &se_cfg, config_t &cfg, int energy);

#endif // __CONFIG_H__