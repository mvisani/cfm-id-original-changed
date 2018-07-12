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

#include "EmModel.h"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include "ModelBase.h"

void ModelBase::initComms() {
    // Initialise the communicator
    int mpi_rank, mpi_nump;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);
    if (mpi_rank == MASTER)
        comm = new MasterComms();
    else
        comm = new WorkerComms();
}

void ModelBase::writeStatus(const char *msg) {

    std::ofstream out;
    out.open(status_filename.c_str(), std::ios_base::out | std::ios_base::app);
    out << msg << std::endl;
    out.close();
}

void ModelBase::writeParamsToFile(std::string &filename) {
    param->saveToFile(filename);
}

void ModelBase::zeroUnusedParams() {

    unsigned int i;
    for (i = 0; i < param->getNumWeights(); i++) {
        if (((MasterComms *) comm)->master_used_idxs.find(i) ==
            ((MasterComms *) comm)->master_used_idxs.end())
            param->setWeightAtIdx(0.0, i);
    }
}

void ModelBase::getEnergiesLevels(std::vector<unsigned int> &energies) {
    unsigned int energy;
    int prev_energy = -1;
    for (unsigned int d = 0; d < cfg->model_depth; d++) {
        energy = cfg->map_d_to_energy[d];
        if (energy != prev_energy)
            energies.push_back(energy);
        prev_energy = energy;
    }
}

void ModelBase::setMiniBatchFlags(std::vector<int> &minibatch_flags, int num_batch) {

    int idx = 0;
    for(auto & flag: minibatch_flags){
        flag = idx;
        idx ++;
        if(idx == num_batch)
            idx = 0;
    }
    std::random_shuffle(minibatch_flags.begin(),minibatch_flags.end());
}

Solver *ModelBase::getSolver(int ga_method, double learning_rate) const {
    Solver *solver;
    switch (ga_method) {
        case USE_ADAM_FOR_GA:
            solver = new Adam(param->getNumWeights(),
                              learning_rate,
                              cfg->ga_adam_beta_1,
                              cfg->ga_adam_beta_2,
                              cfg->ga_adam_eps);

            break;
        case USE_ADAMW_FOR_GA:
            solver = new AdamW(param->getNumWeights(),
                              learning_rate,
                              cfg->ga_adam_beta_1,
                              cfg->ga_adam_beta_2,
                              cfg->ga_adam_eps,
                              cfg->ga_adamw_w);
            break;
        case USE_AMSGRAD_FOR_GA:
            solver = new AMSgrad(param->getNumWeights(),
                                 learning_rate,
                                 cfg->ga_adam_beta_1,
                                 cfg->ga_adam_beta_2,
                                 cfg->ga_adam_eps);

            break;
        case USE_ADADELTA_FOR_GA:
            solver = new Adadelta(param->getNumWeights(),
                                  learning_rate,
                                  cfg->ga_adadelta_rho,
                                  cfg->ga_adam_eps);
            break;
        case USE_MOMENTUM_FOR_GA:
            solver = new Momentum(param->getNumWeights(),
                                  learning_rate,
                                  cfg->ga_momentum);
            break;
        case USE_SGD_FOR_GA:
        default:
            solver = new Sgd(learning_rate);

    }
    return solver;
}

void ModelBase::initParams() {
    switch (cfg->param_init_type){
        case PARAM_FULL_ZERO_INIT:
            param->fullZeroInit();
            break;
        case PARAM_ZERO_INIT:
            param->zeroInit();
            break;
        case PARAM_NORMAL_INIT:
            param->randomNormalInit();
            break;
        case PARAM_RANDOM_INIT:
        default:
            param->randomUniformInit();
    }
}