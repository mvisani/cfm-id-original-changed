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

#include <boost/container/vector.hpp>
#include <string>
#include <vector>

// Default Values
static const std::string APP_STRING = "CFM-ID";

// Fragmentation depth in the initial fragmentation graph
static const int DEFAULT_FRAGGRAPH_DEPTH = 2;

// Fragmentation depth in the model
static const int DEFAULT_MODEL_DEPTH = 2;

// Observation variance (parameters to the peak gaussian distribution -
//  indicative of mass spec mass accuracy)
static const double DEFAULT_ABS_MASS_TOL = 0.01;
static const double DEFAULT_PPM_MASS_TOL = 10.0;

// Threshold below which the gradient needs to drop to end gradient ascent (in UpdateParameters)
static const double DEFAULT_GA_CONVERGE_THRESH = 0.0001;

// Regularization Constant
static const double DEFAULT_LAMBDA = 0.01;

// GA methods
static const int USE_SGD_FOR_GA       = 0;
static const int USE_MOMENTUM_FOR_GA  = 1;
static const int USE_ADAM_FOR_GA      = 2;
static const int USE_ADAMW_FOR_GA     = 3;
static const int USE_ADADELTA_FOR_GA  = 4;
static const int USE_ADABELIEF_FOR_GA = 5;

// Use ADAM to do gradient ascent
static const int DEFAULT_USE_ADAM_FOR_GA = USE_ADAM_FOR_GA;

// LR Decay methods
static const int USE_NO_DECAY      = 0;
static const int USE_DEFAULT_DECAY = 1;
static const int USE_EXP_DECAY     = 2;
static const int USE_STEP_DECAY    = 3;

// Decay
static const float DEFAULT_DECAY_RATE           = 0.1;
static const float DEFAULT_EXP_DECAY_K          = 0.001;
static const float DEFAULT_STEP_DECAY_DROP      = 0.5;
static const int DEFAULT_STEP_DECAY_EPOCHS_DROP = 10;

// Learning Rate
static const float DEFAULT_LEARNING_RATE = 0.001;

// Momentum term used in Gradient Ascent
static const float DEFAULT_MOMENTUM_ALPHA = 0.1;

// ADAM
static const float DEFAULT_ADAM_BETA_1 = 0.9;
static const float DEFAULT_ADAM_BETA_2 = 0.999;

// AdaDelta
static const float DEFAULT_ADADELTA_LEARNING_RATE = 1.0;
static const float DEFAULT_ADADELTA_RHO           = 0.95;

// EPS
static const float DEFAULT_EPS = 1e-8;

// Max by which Q can change between iterations to call convergence
static const float DEFAULT_EM_CONVERGE_THRESH = 0.001;

// Number of times to do a random restart of EM
static const int DEFAULT_NUM_EM_RESTARTS = 3;

static const int POSITIVE_ESI_IONIZATION_MODE = 1;
static const int NEGATIVE_ESI_IONIZATION_MODE = 2;
static const int POSITIVE_EI_IONIZATION_MODE  = 3;

static const int DEFAULT_IONIZATION_MODE   = POSITIVE_ESI_IONIZATION_MODE;
static const int DEFAULT_INCLUDE_ISOTOPES  = 0;
static const double DEFAULT_ISOTOPE_THRESH = 1.0; // On a normalized scale where the highest peak has height 100.0
static const int DEFAULT_INCLUDE_H_LOSSES  = 0;   // For backwards compatibility with old ESI-MS/MS models
static const int DEFAULT_INCLUDE_PRECURSOR_H_LOSSES_ONLY = 0;

static const int PARAM_FULL_ZERO_INIT      = 1;
static const int PARAM_ZERO_INIT           = 2;
static const int PARAM_RANDOM_INIT         = 3;
static const int PARAM_NORMAL_INIT         = 4;
static const int NN_PARAM_VAR_SCALING_INIT = 5;

static const int PARAM_DEFAULT_INIT = PARAM_RANDOM_INIT;

static const int DEFAULT_ALLOW_FRAG_DETOURS = 1;
static const int DEFAULT_MAX_RING_BREAKS    = 2;

static const int MAX_BREAKABLE_RING_SIZE = 8;

// Mode for writing spectra to output
static const int NO_OUTPUT_MODE         = 0;
static const int MSP_OUTPUT_MODE        = 1;
static const int MGF_OUTPUT_MODE        = 2;
static const int SINGLE_TXT_OUTPUT_MODE = 3;

// Maximum additional electron pairs to either side during fragmentation
static const int MAX_E_MOVE = 4;

static const int LINEAR_THETA_FUNCTION     = 1;
static const int NEURAL_NET_THETA_FUNCTION = 2;
static const int DEFAULT_THETA_FUNCTION    = LINEAR_THETA_FUNCTION;

static const int LINEAR_NN_ACTIVATION_FUNCTION            = 0;
static const int RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION = 1;
static const int RELU_NN_ACTIVATION_FUNCTION              = 2;
static const int LEAKY_RELU_NN_ACTIVATION_FUNCTION        = 3;

static const int DEFAULT_NN_ACTIVATION_FUNCTION = RELU_AND_NEG_RLEU_NN_ACTIVATION_FUNCTION;

static const double DEFAULT_GA_MOMENTUM = 0.9;

static const int DEFAULT_GA_MINIBATCH_NTH_SIZE = 1;
static const int DEFAULT_GA_MAX_ITERATIONS     = 20;
static const int DEFAULT_EM_MAX_ITERATIONS     = 100;

static const int NORMAL_OBS_FUNCTION  = 1;
static const int UNIFORM_OBS_FUNCTION = 2;
static const int DEFAULT_OBS_FUNCTION = NORMAL_OBS_FUNCTION;

// Timeout settings
static const int DEFAULT_FRAGGRAPH_COMPUTE_TIMEOUT_IN_SECS = -1; //-1 is no timeout.

// Random Sample settings
static const int DEFAULT_USE_BEST_Q_IN_GA       = 0;
static const int USE_NO_SAMPLING                = 0;
// Random selection
static const int USE_RANDOM_SAMPLING            = 1;
// Random  walk selection
static const int USE_GRAPH_RANDOM_WALK_SAMPLING = 2;
// Selection child node base on intensity difference
static const int USE_DIFFERENCE_SAMPLING_BFS_CO = 3;
// Selection all node (full path) base on intensity difference
static const int USE_DIFFERENCE_SAMPLING_BFS    = 4;

// Configuration
struct config_t {
	// unsigned int num_threads;

	// Fragment Graph Configuration
	int fg_depth;
	int allow_frag_detours;
	int max_ring_breaks;
	int include_h_losses;
	int include_precursor_h_losses_only;

	// Use Single Energy CFM (rather than Combined Energy)
	int ionization_mode;
	int include_isotopes;
	double isotope_thresh;
	std::string isotope_pattern_file;

	// For sampling
	bool ga_reset_sampling;
	int ga_sampling_method;
	int ga_sampling_method2;
	int ga_sampling_max_selection;
	int ga_diff_sampling_peak_num;
	double ga_diff_sampling_difference;

	// Model Level Configuration
	unsigned int model_depth; // Total Depth
	std::vector<int> spectrum_depths;
	std::vector<double> spectrum_weights;
	double abs_mass_tol;
	double ppm_mass_tol;
	std::vector<int> map_d_to_energy;     // Derived parameter
	std::vector<int> dv_spectrum_depths;  // Either a direct copy, or interpolated values.
	std::vector<int> dv_spectrum_indexes; // Index of each spectrum in the list of spectra in MolData
	std::vector<double> dv_spectrum_weights;
	int obs_function; // Function used for the observation function

	// Theta function Configuration (Linear or Neural Net)
	int theta_function;
	// Neural Net Configuration
	std::vector<int> theta_nn_hlayer_num_nodes;
	std::vector<int> theta_nn_layer_act_func_ids;
	std::vector<float> nn_layer_dropout_probs;
	boost::container::vector<bool> nn_layer_is_frozen_flags;

	// EM Configuration
	double em_converge_thresh;
	int num_em_restarts;
	int use_lower_energy_params_for_init;
	int param_init_type;

	// Gradient Ascent Configuration
	double lambda; // Regularization n
	int ga_method;
	int ga_decay_method;
	double ga_converge_thresh;

	float starting_step_size;
	float ending_step_size;
	int ga_max_iterations;
	int ga_use_best_q;
	int em_max_iterations;

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
	bool ga_adam_use_amsgrad = false;
	// For ADAMW
	double ga_adamw_w;

	// For adadelta
	double ga_adadelta_rho;

	int update_bias_first;
	int ga_minibatch_nth_size;

	int fragraph_compute_timeout_in_secs;

	bool disable_cross_val_metrics;
	bool disable_training_metrics;
	bool disable_cpu_usage_metrics;

	int em_no_progress_count;
	int ga_no_progress_count;

	bool collected_all_used_idx;
	bool allow_intermediate_peak;
	bool allow_cyclization;

	bool use_log_scale_peak;
	bool use_iterative_fg_gen;

	// default post-processing settings
	int default_predicted_peak_min;
	int default_predicted_peak_max;
	int default_mz_decimal_place;
	double default_predicted_min_intensity;
	double default_postprocessing_energy;

	double intensity_msg_weight;
};

// argv_zero: argv[0], use for extract binary folder
void initDefaultConfig(config_t &cfg, char *argv_zero = nullptr);

// filename: config filename
// argv_zero: argv[0], use for extract binary folder
void initConfig(config_t &cfg, std::string &filename, char *argv_zero, bool report_all);

void initDerivedConfig(config_t &cfg, int energy = -1);

void initSingleEnergyConfig(config_t &se_cfg, config_t &cfg, int energy);

void printSamplingConfig(const int &sampling_method, config_t &cfg);

#endif // __CONFIG_H__
