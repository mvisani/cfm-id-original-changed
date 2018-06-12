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
    //int num_per_batch = minibatch_flags.size() / num_batch;

    // The flags are initialized to 1's for selectable molecules, so set unwanted
    // molecule flags to 0
    /*int num_mols = initialized_minibatch_flags.size();
    std::vector<int> idxs(num_mols);
    int count = 0;
    for (int i = 0; i < num_mols; i++)
        if (initialized_minibatch_flags[i])
            idxs[count++] = i;
    idxs.resize(count);


    shuffle(idxs.begin(), idxs.end(), util_rng);

    int num_minibatch_mols =
            (num_mols + cfg->ga_minibatch_nth_size - 1) / cfg->ga_minibatch_nth_size;
    for (int i = num_minibatch_mols; i < idxs.size(); i++)
        initialized_minibatch_flags[idxs[i]] = 0;*/
}

Solver *ModelBase::getSolver(int ga_method, double learning_rate) const {
    Solver *solver;
    switch (ga_method) {
        case USE_ADAM_FOR_GA:
            solver = new Aadm(param->getNumWeights(),
                              learning_rate,
                              cfg->ga_adam_beta_1,
                              cfg->ga_adam_beta_2,
                              cfg->ga_eps);

            break;
        case USE_AMSGRAD_FOR_GA:
            solver = new AMSgrad(param->getNumWeights(),
                                 learning_rate,
                                 cfg->ga_adam_beta_1,
                                 cfg->ga_adam_beta_2,
                                 cfg->ga_eps);

            break;
        case USE_ADADELTA_FOR_GA:
            solver = new Adadelta(param->getNumWeights(),
                                  learning_rate,
                                  cfg->ga_adadelta_rho,
                                  cfg->ga_eps);
            break;
        case USE_MOMENTUM_FOR_GA:
        default:
            solver = new Momentum(param->getNumWeights(),
                                  learning_rate,
                                  cfg->ga_momentum);
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