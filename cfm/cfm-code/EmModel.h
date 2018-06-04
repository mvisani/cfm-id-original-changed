/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.h
#
# Description: 	Class to apply Expectation Maximization algorithm to derive 
#				model parameters.
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

#ifndef __EM_TRAIN_H__
#define __EM_TRAIN_H__

static const int MAX_EM_ITERATIONS = 100;

#include "Config.h"
#include "MolData.h"
#include "Param.h"
#include "NNParam.h"
#include "IPFP.h"
#include "Solver.h"
#include "ModelBase.h"

//Sufficient stats, access by transition index for a given molecule
typedef std::vector<double> suft_t;

struct suft_counts_t {
    //Access each by molecule index
    std::vector<suft_t> values;
};

class EmComputationException : public std::exception {

    virtual const char *what() const noexcept override {
        return "Error during EM, unable to proceed.";
    }
};


class EmModel : public ModelBase {
public:
    //Constructor
    //Note: To include the group in the status filename, include _GRP_ in the name
    EmModel(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename,
       std::string initial_params_filename = "");

    ~EmModel();

    //Run the EM algorithm on the supplied data (except the specified group),
    //return the final likelihood value.
    double
    trainModel(std::vector<MolData> &molDataSet, int group, std::string &out_param_filename, int energy_level) override;

    //This is public so the test can access it....there must be a better way?
    virtual double computeAndAccumulateGradient(double *grads, int molidx, MolData &moldata, suft_counts_t &suft,
                                                bool record_used_idxs_only, std::set<unsigned int> &used_idxs,
                                                int sampling_method = 0);

    virtual double computeQ(int molidx, MolData &moldata, suft_counts_t &suft);

protected:

    //Further virtual functions
    virtual void computeThetas(MolData *moldata);

    //Use to note when unused parameters have been zeroed (so we don't

    //Initialise sufficient statistics
    void initSuft(suft_counts_t &suft, std::vector<MolData> &data);

    //Update sufficient statistics based on beliefs
    void recordSufficientStatistics(suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs);


    //Simple gradient ascent
    double updateParametersGradientAscent(std::vector<MolData> &data, suft_counts_t &suft, double learning_rate,
                                          int sampling_method);

    //Helper functions

    // function to add Regularization term for Q
    // update grads is updated if a ptr is passed
    // nullptr means do not update grads
    virtual double addRegularizersAndUpdateGradient(double *grads);

};

#endif // __EM_TRAIN_H__
