/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.h
#
# Description: 	Class to apply Direact GA algorithm to derive
#				model parameters.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#ifndef CFM_DIRECTMODEL_H
#define CFM_DIRECTMODEL_H

#include "ModelBase.h"

class DirectModel: public ModelBase {

public:
    //This is public so the test can access it....there must be a better way?
    virtual double ComputeAndAccumulateGradient(double *grads, MolData &moldata,
                                                bool record_used_idxs_only, std::set<unsigned int> &used_idxs,
                                                int sampling_method);

    virtual double ComputeQ(MolData &moldata);

    double TrainModel(std::vector<MolData> &molDataSet, int group, std::string &out_param_filename) override;

    // function to add Regularization term for Q
    // update grads is updated if a ptr is passed
    // nullptr means do not update grads
    virtual double addRegularizersAndUpdateGradient(double *grads);

protected:
    //Simple gradient ascent
    double updateParametersGradientAscent(std::vector<MolData> &data, double learning_rate, int sampling_method);

};


#endif //CFM_DIRECTMODEL_H
