/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.cpp
#
# Description: 	Class to apply Machine Learning algorithm to derive
#				model parameters.

#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "ModelBase.h"
#include "EmModel.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

void ModelBase::writeStatus(const char *msg) {

	std::ofstream out;
	out.open(status_filename.c_str(), std::ios_base::out | std::ios_base::app);
	out << msg << std::endl;
	out.close();
}

void ModelBase::writeParamsToFile(std::string &filename) { param->saveToFile(filename); }

void ModelBase::getEnergiesLevels(std::vector<unsigned int> &energies) {
	unsigned int energy;
	int prev_energy = -1;
	for (unsigned int d = 0; d < cfg->model_depth; d++) {
		energy = cfg->map_d_to_energy[d];
		if (energy != prev_energy) energies.push_back(energy);
		prev_energy = energy;
	}
}

void ModelBase::setMiniBatchFlags(std::vector<int> &minibatch_flags, int num_batch) {

	int idx = 0;
	for (auto &flag : minibatch_flags) {
		flag = idx;
		idx++;
		if (idx == num_batch) idx = 0;
	}
	std::shuffle(minibatch_flags.begin(), minibatch_flags.end(), util_rng);
}

Solver *ModelBase::getSolver(int ga_method, double learning_rate) const {
	Solver *solver;
	switch (ga_method) {
	case USE_ADAM_FOR_GA:
		solver = new Adam(param->getNumWeights(), learning_rate, cfg->ga_adam_beta_1, cfg->ga_adam_beta_2,
		                  cfg->ga_adam_eps, cfg->ga_adam_use_amsgrad);

		break;
	case USE_ADAMW_FOR_GA:
		solver = new AdamW(param->getNumWeights(), learning_rate, cfg->ga_adam_beta_1, cfg->ga_adam_beta_2,
		                   cfg->ga_adam_eps, cfg->ga_adam_use_amsgrad, cfg->ga_adamw_w);
		break;
	case USE_ADABELIEF_FOR_GA:
		solver = new AdaBelief(param->getNumWeights(), learning_rate, cfg->ga_adam_beta_1, cfg->ga_adam_beta_2,
		                       cfg->ga_adam_eps, cfg->ga_adam_use_amsgrad);

		break;
	case USE_ADADELTA_FOR_GA:
		solver = new AdaDelta(param->getNumWeights(), learning_rate, cfg->ga_adadelta_rho, cfg->ga_adam_eps);
		break;
	case USE_MOMENTUM_FOR_GA: solver = new Momentum(param->getNumWeights(), learning_rate, cfg->ga_momentum); break;
	case USE_SGD_FOR_GA:
	default: solver = new Sgd(learning_rate);
	}
	return solver;
}
