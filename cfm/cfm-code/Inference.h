/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# inference.h
#
# Description: 	Custom implementation of jtree based exact inference for
#				this network (for efficiency).
#				
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __INFER_H__
#define __INFER_H__

#include "MolData.h"
#include "Message.h"
#include "Config.h"

class InvalidInferenceException : public std::exception {

    virtual const char *what() const noexcept {
        return "Error during inference calculation, unable to proceed.";
    }
};

typedef std::vector<std::vector<double> > spec_factor_t;

struct factor_probs_t {
    std::vector<std::vector<double> > tn;    //Transition:  Indexed by transition, then depth
    std::vector<std::vector<double> > ps;    //Persistence: Indexed by fragment, then depth
};

class Inference {

public:
    Inference(const MolData *data, config_t *cfg) : moldata(data), config(cfg) {
        if(data->getFGHeight() <= config->model_depth)
            mol_depth = config->model_depth;
        else
            mol_depth = data->getFGHeight();

        intensity_msg_weight = cfg->intensity_msg_weight;
    };

    void calculateBeliefs(beliefs_t &beliefs, int current_energy);    //Note: Only supports single energy mode
    void runInferenceDownwardPass(std::vector<Message> &down_msgs, int to_depth, int energy);

    void createSpectrumMessage(Message &msg, int energy, Message &down_msg);

private:
    const MolData *moldata;
    config_t *config;
    std::vector<Message> spec_messages;
    int mol_depth;
    double intensity_msg_weight;

    void initTmpFactorProbSizes(factor_probs_t &tmp_log_probs, unsigned int num_frag, unsigned int num_trans,
                                unsigned int model_depth);

    void runInferenceUpwardPass(std::vector<Message> &up_msgs, Message &spec_msg, int current_energy);

    void createMessage(factor_probs_t &tmp_log_probs, Message &m, Message &prev_m, int direction, int depth);

    void passMessage(factor_probs_t &tmp_log_probs, int direction, int depth, Message &m, int energy);

    void
    combineMessagesToComputeBeliefs(beliefs_t &beliefs, std::vector<Message> &down_msgs,
                                    std::vector<Message> &up_msgs, int current_energy);

    void createSpectrumMessageWithIsotopes(Message &msg, int energy, Message &down_msg);

};

#endif // __INFER_H__