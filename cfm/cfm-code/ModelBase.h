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
#ifndef CFM_MODELBASE_H
#define CFM_MODELBASE_H

#include "Comms.h"
#include "Config.h"
#include "IPFP.h"
#include "MolData.h"
#include "Solver.h"

#include "mpi.h"

class ModelBase {

protected:
    //do it more times than we need to)
    void initComms();

    void writeStatus(const char *msg);

    //The feature calculator to use - preconfigured with feature spec
    FeatureCalculator *fc;
    //The current parameters
    boost::shared_ptr<Param> param;
    //Configuration data
    config_t *cfg;
    //Communicator (for exchanging data between threads)
    Comms *comm;
    int unused_zeroed;
    //For writing status messages to a log file
    std::string status_filename;
    int validation_group;
    bool sparse_params;

    void zeroUnusedParams();

    // function to get Eng levels
    void getEnergiesLevels(std::vector<unsigned int> &energies);

    Solver *getSolver(int ga_method, double learning_rate) const;

    virtual void initParams();

    bool initial_params_provided;
public:
    //After running Model, the final params can be written out to file
    virtual void writeParamsToFile(std::string &filename);

    // abstract function
    virtual double
    trainModel(std::vector<MolData> &molDataSet, int group, std::string &out_param_filename, int energy_level) = 0;

    //Select mini batch (exposed publicly for testing...)
    void setMiniBatchFlags(std::vector<int> &minibatch_flags, int num_batch);
};


#endif //CFM_MODELBASE_H
