/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# inference.cpp
#
# Description: 	Custom implementation of jtree based exact inference for
#				this exact network (for efficiency).
#				
#
# Copyright (c) 2013, Felicity Allen ,2017 ,Fei Wang
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Inference.h"
#include "Util.h"

#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions.hpp>

//Big log negative value that will effectively result in zero probability
static const double NULL_PROB = -A_BIG_DBL;
static const int DOWN = 0;
static const int UP = 1;

void Inference::calculateBeliefs(beliefs_t &beliefs, int current_energy) {
    /*
     * Disabled beacsuse we only havs se-cfm now
    if (!config->use_single_energy_cfm) {
        std::cout << "Error: Use of inference calculateBeliefs only valid for single energy model. Try IPFP instead."
                  << std::endl;
        throw InvalidInferenceException();
    }*/

    //Pass the messages down
    std::vector<Message> down_msgs, up_msgs;
    runInferenceDownwardPass(down_msgs, mol_depth, current_energy);

    //Create Spectrum Message
    Message spec_msg;
    createSpectrumMessage(spec_msg, current_energy, down_msgs[mol_depth - 1]);

    //Apply IPFP modification to message to account for marginal observation
    Message modified_msg;
    modified_msg.reset(spec_msg.size());
    Message::const_iterator it = spec_msg.begin();
    for (; it != spec_msg.end(); ++it) {
        int i = it.index();
        modified_msg.addToIdx(i, spec_msg.getIdx(i) - down_msgs[mol_depth - 1].getIdx(i));
    }

    //Pass the messages back up
    runInferenceUpwardPass(up_msgs, modified_msg, current_energy);

    //Combine the messages and initial probabilities to compute the beliefs
    combineMessagesToComputeBeliefs(beliefs, down_msgs, up_msgs, current_energy);
}


void Inference::initTmpFactorProbSizes(factor_probs_t &tmp_log_probs, unsigned int num_frag, unsigned int num_trans,
                                       unsigned int model_depth) {

    //Initialise persistence vector sizes
    tmp_log_probs.ps.resize(num_frag);
    for (unsigned int i = 0; i < num_frag; i++)
        tmp_log_probs.ps[i].resize(model_depth - 1);

    //Initialise transition vector sizes
    tmp_log_probs.tn.resize(num_trans);
    for (unsigned int i = 0; i < num_trans; i++)
        tmp_log_probs.tn[i].resize(model_depth - 1);
}


void Inference::runInferenceDownwardPass(std::vector<Message> &down_msgs, int to_depth, int energy) {

    //Initialise the messages
    down_msgs.resize(mol_depth);

    //Create the tmp factor probs
    factor_probs_t tmp_log_probs;
    initTmpFactorProbSizes(tmp_log_probs, moldata->getNumFragments(), moldata->getNumTransitions(), mol_depth);

    //Factor (F0,F1) => Create F1 Message
    const std::vector<int> *tmap = &((*moldata->getFromIdTMap())[0]);
    std::vector<int>::const_iterator it = tmap->begin();
    down_msgs[0].reset(moldata->getNumFragments());
    down_msgs[0].addToIdx(0, moldata->getLogPersistenceProbForIdx(energy, 0));

    for (; it != tmap->end(); ++it) {
        auto t = moldata->getTransitionAtIdx(*it);
        down_msgs[0].addToIdx(t->getToId(), moldata->getLogTransitionProbForIdx(energy, *it));
    }

    //Update Factor (F1,F2) => Create Message F2 => Update Factor (F2,F3) ...etc as per MODEL_DEPTH
    for (int i = 0; i < to_depth - 1; i++) {
        passMessage(tmp_log_probs, DOWN, i, down_msgs[i], energy);
        createMessage(tmp_log_probs, down_msgs[i + 1], down_msgs[i], DOWN, i);
    }
}

void Inference::runInferenceUpwardPass(std::vector<Message> &up_msgs, Message &spec_msg, int current_energy) {


    //Initialise the messages
    up_msgs.resize(mol_depth);

    //Create the tmp factor probs
    factor_probs_t tmp_log_probs;
    initTmpFactorProbSizes(tmp_log_probs, moldata->getNumFragments(), moldata->getNumTransitions(), mol_depth);

    //Apply spectrum message
    up_msgs[mol_depth - 1] = spec_msg;

    //Update Factor (Fd-1,Fd) => Create Message Fd-1  ...etc as per MODEL_DEPTH
    for (int i = mol_depth - 2; i >= 0; i--) {
        // make sure we are not go to far
        //(i+1 <= current_energy) ? config->map_d_to_energy[i + 1] : current_energy;
        passMessage(tmp_log_probs, UP, i, up_msgs[i + 1], current_energy);
        createMessage(tmp_log_probs, up_msgs[i], up_msgs[i + 1], UP, i);
    }
}

void Inference::createMessage(factor_probs_t &tmp_log_probs, Message &m, Message &prev_m, int direction, int depth) {

    m.reset(moldata->getNumFragments());
    for (unsigned int id = 0; id < moldata->getNumFragments(); id++) {
        //Going up or down?
        const std::vector<int> *tmap;
        if (direction == DOWN)
            tmap = &((*moldata->getToIdTMap())[id]);
        else
            tmap = &((*moldata->getFromIdTMap())[id]);

        //Marginalize out the upper/lower variable
        //the log_sum of ps term means the chance to remain as this fragment at this step?
        double log_sum = -A_BIG_DBL;
        if (prev_m.getIdx(id) > -A_BIG_DBL)
            log_sum = tmp_log_probs.ps[id][depth];

        for (auto itt = tmap->begin(); itt != tmap->end(); ++itt) {
            auto t = moldata->getTransitionAtIdx(*itt);
            // if perm[fragment_id] is not super small
            // nothing should happen
            // else, prob to next depth - ps term at this depth  * trans term at this depth
            if ((direction == DOWN && prev_m.getIdx(t->getFromId()) > -A_BIG_DBL) ||
                (direction == UP && prev_m.getIdx(t->getToId()) > -A_BIG_DBL)) {
                log_sum = logAdd(log_sum, tmp_log_probs.tn[*itt][depth]);
            }
        }
        if (log_sum > -A_BIG_DBL)
            m.addToIdx(id, log_sum);
    }
}

void Inference::passMessage(factor_probs_t &tmp_log_probs, int direction, int depth, Message &m, int energy) {

    Message::const_iterator it = m.begin();
    for (; it != m.end(); ++it) {

        unsigned int idx = it.index();

        //Apply to the persistence term (idx -> idx)
        tmp_log_probs.ps[idx][depth] = moldata->getLogPersistenceProbForIdx(energy, idx) + m.getIdx(idx);

        //Apply to all the other transitions applicable for this message element
        const std::vector<int> *tmap;
        if (direction == DOWN)
            tmap = &((*moldata->getFromIdTMap())[idx]);
        else
            tmap = &((*moldata->getToIdTMap())[idx]);
        std::vector<int>::const_iterator itt = tmap->begin();
        for (; itt != tmap->end(); ++itt)
            tmp_log_probs.tn[*itt][depth] = moldata->getLogTransitionProbForIdx(energy, *itt) + m.getIdx(idx);
    }
}

void Inference::combineMessagesToComputeBeliefs(beliefs_t &beliefs, std::vector<Message> &down_msgs,
                                                std::vector<Message> &up_msgs, int current_energy) {

    std::vector<double> norms(mol_depth);

    //Compute Persistence Beliefs (and track norms)
    beliefs.ps.resize(moldata->getNumFragments());
    for (unsigned int i = 0; i < moldata->getNumFragments(); i++) {
        beliefs.ps[i].resize(mol_depth);

        for (unsigned int d = 0; d < mol_depth; d++) {

            double tmp;
            if ((d == 0 && i == 0) || (d > 0 && down_msgs[d - 1].getIdx(i) > -A_BIG_DBL)) {
                tmp = moldata->getLogPersistenceProbForIdx(current_energy, i);
                tmp += up_msgs[d].getIdx(i);
                if (d > 0) tmp += down_msgs[d - 1].getIdx(i);
                if (i == 0) norms[d] = tmp;
                else norms[d] = logAdd(norms[d], tmp);
            } else tmp = NULL_PROB;
            beliefs.ps[i][d] = tmp;
        }
    }

    //Compute Transition Beliefs (and track norms)
    beliefs.tn.resize(moldata->getNumTransitions());
    for (unsigned int i = 0; i < moldata->getNumTransitions(); i++) {
        auto t = moldata->getTransitionAtIdx(i);
        beliefs.tn[i].resize(mol_depth);

        for (unsigned int d = 0; d < mol_depth; d++) {
            double tmp;
            if ((d == 0 && t->getFromId() == 0) || (d > 0 && down_msgs[d - 1].getIdx(t->getFromId()) > -A_BIG_DBL)) {
                tmp = moldata->getLogTransitionProbForIdx(current_energy, i);
                // since we are using log, log(prob a * prob b) = log(prob a) + log(prob b)
                tmp += up_msgs[d].getIdx(t->getToId());
                if (d > 0)
                    tmp += down_msgs[d - 1].getIdx(t->getFromId());
                norms[d] = logAdd(norms[d], tmp);
            } else tmp = NULL_PROB;
            beliefs.tn[i][d] = tmp;
        }
    }

    //Normalise
    for (unsigned int i = 0; i < beliefs.tn.size(); i++) {
        for (unsigned int d = 0; d < mol_depth; d++) {
            beliefs.tn[i][d] -= norms[d];
        }
    }
    for (unsigned int i = 0; i < beliefs.ps.size(); i++) {
        for (unsigned int d = 0; d < mol_depth; d++) {
            beliefs.ps[i][d] -= norms[d];
        }
    }
}

void Inference::createSpectrumMessage(Message &msg, int energy, Message &down_msg) {

    const Spectrum *spectrum = moldata->getSpectrum(energy);

    if (moldata->hasIsotopesIncluded())
        createSpectrumMessageWithIsotopes(msg, energy, down_msg);
    else {
        //Store normpdf( pk mass, ion mass, sigma*sqrt2 )
        static const double pi = boost::math::constants::pi<double>();
        int num_fragments = moldata->getNumFragments();
        msg.reset(num_fragments);

        Spectrum::const_iterator pk = spectrum->begin();
        for (; pk != spectrum->end(); ++pk) {
            //This will store the fragment probabilites for this single peak
            //(each message needs to be normalised independently before combining)
            Message peak_msg;
            peak_msg.reset(num_fragments);

            //Set peak sigma based on the specified mass tolerances and the peak of interest
            //(assume tolerances cut things off at approx 3 x std deviation)
            double peak_sigma = 0.33 * getMassTol(config->abs_mass_tol, config->ppm_mass_tol, pk->mass);
            double norm = -0.5 * std::log(4 * pi * peak_sigma * peak_sigma);
            double denom = 0.25 / (peak_sigma * peak_sigma);
            for (unsigned int j = 0; j < num_fragments; j++) {

                //Disallow fragments we can't get to in the model.
                if (down_msg.getIdx(j) < -10000.0)
                    continue;

                //Disallow fragments from explaining distant peaks
                //(problematic in the absence of a better explanation)
                const Fragment *fgt = moldata->getFragmentAtIdx(j);
                double mass_diff = fabs(fgt->getMass() - pk->mass);
                if (mass_diff > 3 * peak_sigma)
                    continue;

                double sq_mass_diff =  mass_diff * mass_diff;
                if (config->obs_function == UNIFORM_OBS_FUNCTION)
                    sq_mass_diff  =  0.0;

                //The down message is applied here to weight competing
                //fragments for the same peak based on current probability
                //estimates.
                peak_msg.addToIdx(j, norm - denom * sq_mass_diff + down_msg.getIdx(j));
            }
            msg.addWeightedMessage(peak_msg, intensity_msg_weight * pk->intensity);
        }
    }
}

void Inference::createSpectrumMessageWithIsotopes(Message &msg, int energy, Message &down_msg) {

    const Spectrum *spectrum = moldata->getSpectrum(energy);

    //Store normpdf( pk mass, ion mass, sigma*sqrt2 )
    static const double pi = boost::math::constants::pi<double>();
    int num_fragments = moldata->getNumFragments();
    msg.reset(num_fragments);

    Spectrum::const_iterator pk = spectrum->begin();
    for (; pk != spectrum->end(); ++pk) {

        Message peak_msg;    //This will store the fragment probabilites for this single peak
        //(each message needs to be normalised independently before combining)
        peak_msg.reset(num_fragments);

        //Set peak sigma based on the specified mass tolerances and the peak of interest
        //(assume tolerances cut things off at approx 3 x std deviation)
        double peak_sigma = 0.33 * getMassTol(config->abs_mass_tol, config->ppm_mass_tol, pk->mass);
        double norm = -0.5 * std::log(4 * pi * peak_sigma * peak_sigma);
        double denom = 0.25 / (peak_sigma * peak_sigma);
        for (unsigned int j = 0; j < num_fragments; j++) {
            if (down_msg.getIdx(j) < -10000.0) continue;    //Disallow fragments we can't get to in the model.
            const Fragment *fgt = moldata->getFragmentAtIdx(j);

            const Spectrum *iso_spectrum = fgt->getIsotopeSpectrum();
            Spectrum::const_iterator ipk = iso_spectrum->begin();
            for (; ipk != iso_spectrum->end(); ++ipk) {
                double mass_diff = fabs(ipk->mass - pk->mass);
                if (mass_diff > 3 * peak_sigma) continue;    //Disallow fragments from explaining distant peaks
                //(problematic in the absence of a better explanation)

                double sq_mass_diff =  mass_diff * mass_diff;
                if (config->obs_function == UNIFORM_OBS_FUNCTION)
                    sq_mass_diff  =  0.0;

                peak_msg.addToIdx(j,
                                  norm - denom * sq_mass_diff + std::log(intensity_msg_weight* ipk->intensity) + down_msg.getIdx(j));
                //The down message is applied here to weight competing
                //fragments for the same peak based on current probability
                //estimates.
            }
        }
        msg.addWeightedMessage(peak_msg, intensity_msg_weight * pk->intensity);
        //msg.addWeightedMessage(peak_msg, 0.01* std::log(1.0 + pk->intensity));
    }

}