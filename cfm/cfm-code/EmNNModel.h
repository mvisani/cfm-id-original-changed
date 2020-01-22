/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.h
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

#ifndef __EM_NN_TRAIN_H__
#define __EM_NN_TRAIN_H__

#include "EmModel.h"
#include "NNParam.h"

class EmNNModel : public EmModel {
public:
    //Constructor
    //Note: To include the group in the status filename, include _GRP_ in the name
    EmNNModel(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename,
              std::string initial_params_filename = "");

    //This is public so the test can access it....there must be a better way?
    int computeAndAccumulateGradient(float *grads, int mol_idx, MolData &mol_data, suft_counts_t &suft,
                                     int sampling_method, unsigned int energy) override;

    void collectUsedIdx(MolData &mol_data, std::set<unsigned int> &used_idxs, unsigned int energy) override;

    double computeLogLikelihoodLoss(int molidx, MolData &moldata, suft_counts_t &suft, unsigned int energy) override;

    double getRegularizationTerm(unsigned int energy) override;

    void updateGradientForRegularizationTerm(float *grads, unsigned int energy) override;

    void writeParamsToFile(std::string &filename) override;

private:
    //The current parameters
    boost::shared_ptr<NNParam> nn_param;

    void computeThetas(MolData *moldata) override;
};

#endif // __EM_NN_TRAIN_H__
