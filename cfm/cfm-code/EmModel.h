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
#include "Solver.h"
#include "ModelBase.h"

//Sufficient stats, access by transition index for a given molecule
typedef std::vector<double> suft_t;

struct suft_counts_t {
    //Access each by molecule index
    std::vector<suft_t> values;
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
    virtual int computeAndAccumulateGradient(float *grads, int mol_idx, MolData &mol_data, suft_counts_t &suft,
                                             int sampling_method, unsigned int energy);

    virtual void collectUsedIdx(MolData &mol_data, std::set<unsigned int> &used_idxs, unsigned int energy);

    virtual double computeLogLikelihoodLoss(int molidx, MolData &moldata, suft_counts_t &suft, unsigned int energy);

protected:

    //Further virtual functions
    virtual void computeThetas(MolData *moldata);

    //Use to note when unused parameters have been zeroed (so we don't

    //Initialise sufficient statistics
    void initSuft(suft_counts_t &suft, std::vector<MolData> &data);

    //Update sufficient statistics based on beliefs
    void recordSufficientStatistics(suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs,
                                    unsigned int energy);


    //Simple gradient ascent
    double updateParametersGradientAscent(std::vector<MolData> &data, suft_counts_t &suft, double learning_rate,
                                          int sampling_method, unsigned int energy);

    double computeLoss(std::vector<MolData> &data, suft_counts_t &suft, unsigned int energy);

    //Helper functions

    // function to add Regularization term for Q
    // update grads is updated if a ptr is passed
    // nullptr means do not update grads
    virtual double getRegularizationTerm(unsigned int energy);

    virtual void updateGradientForRegularizationTerm(float *grads, unsigned int energy);

    void getSubSampledTransitions(MolData &moldata, int sampling_method, unsigned int energy,
                                  std::set<int> &selected_trans_id) const;
    void
    computeMetrics(int energy_level, std::vector<MolData, std::allocator<MolData>>::iterator &itdata,
                       double &jaccard, double &w_jaccard);

    double getUpdatedLearningRate(double learning_rate, int iter) const;

    void updateTrainingParams(double loss, double prev_loss, double loss_ratio, float &learning_rate,
                              int &sampling_method,
                              int &count_no_progress) const;

    bool withinGradOffset(int idx, unsigned energy) {
        unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();
       return (grad_offset * energy) <= idx && (grad_offset * (energy+1) > idx);
    };
};

#endif // __EM_TRAIN_H__
