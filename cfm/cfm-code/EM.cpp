/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.cpp
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

#include "EM.h"
#include "mpi.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

// Init static members
std::random_device   EM::m_rd;
std::mt19937         EM::m_rng(EM::m_rd());
std::uniform_real_distribution<double> EM::m_uniform_dist(0, 1.0);

EM::EM(config_t *a_cfg, FeatureCalculator *an_fc,
       std::string &a_status_filename, std::string initial_params_filename) {
    cfg = a_cfg;
    fc = an_fc;
    status_filename = a_status_filename;
    initComms();
    int num_energies_to_include = cfg->map_d_to_energy.back() + 1;
    if (initial_params_filename == "") {
        param = boost::shared_ptr<Param>(
                new Param(fc->getFeatureNames(), num_energies_to_include));
        initial_params_provided = false;
        comm->printToMasterOnly("EM: No initial params provided");
    } else {
        param = boost::shared_ptr<Param>(new Param(initial_params_filename));
        while (param->getNumEnergyLevels() < num_energies_to_include)
            param->appendRepeatedPrevEnergyParams();
        initial_params_provided = true;
        std::string msg =
                "EM: Initial params provided from " + initial_params_filename;
        comm->printToMasterOnly(msg.c_str());
    }
    sparse_params = true;
}

void EM::initComms() {
    // Initialise the communicator
    int mpi_rank, mpi_nump;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);
    if (mpi_rank == MASTER)
        comm = new MasterComms();
    else
        comm = new WorkerComms();
}

EM::~EM() { delete comm; }

void EM::writeStatus(const char *msg) {

    std::ofstream out;
    out.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
    out << msg << std::endl;
    out.close();
}

void EM::writeParamsToFile(std::string &filename) {
    param->saveToFile(filename);
}

void EM::initParams() {
    if (cfg->em_init_type == PARAM_FULL_ZERO_INIT)
        param->fullZeroInit();
    else if (cfg->em_init_type == PARAM_ZERO_INIT)
        param->zeroInit();
    else
        param->randomInit();
}

void EM::computeThetas(MolData *moldata) {
    moldata->computeTransitionThetas(*param);
}

double EM::run(std::vector<MolData> &data, int group,
               std::string &out_param_filename) {

    unused_zeroed = 0;
    int iter = 0;

    if (!initial_params_provided)
        initParams();
    comm->broadcastInitialParams(param.get());
    validation_group = group;

    // Write the initialised params to file (we may get want to reload and use
    // with saved suft state, even before updating)
    std::string init_out_param_filename = out_param_filename + "_init";
    if (comm->isMaster())
        writeParamsToFile(init_out_param_filename);

    // EM
    iter = 0;
    double Q = 0.0;
    double prevQ = -DBL_MAX;
    double bestQ = -DBL_MAX;

    int count_no_progress = 0;
    while (iter < MAX_EM_ITERATIONS) {

        std::string iter_out_param_filename =
                out_param_filename + "_" + boost::lexical_cast<std::string>(iter);

        std::string msg = "EM Iteration " + boost::lexical_cast<std::string>(iter);
        if (comm->isMaster())
            writeStatus(msg.c_str());
        comm->printToMasterOnly(msg.c_str());

        time_t before, after;
        std::vector<MolData>::iterator itdata;

        // Reset sufficient counts
        suft_counts_t suft;
        initSuft(suft, data);

        int num_converged = 0, num_nonconverged = 0;
        int tot_numc = 0, total_numnonc = 0;
        before = time(nullptr);

        // Do the inference part (E-step)
        itdata = data.begin();
        for (int molidx = 0; itdata != data.end(); ++itdata, molidx++) {

            if (!itdata->hasComputedGraph())
                continue; // If we couldn't compute it's graph for some reason..
            if (itdata->hasEmptySpectrum()) {
                std::cout << "Warning: No peaks with explanatory fragment found for "
                          << itdata->getId() << ", ignoring this input molecule."
                          << std::endl;
                continue; // Ignore any molecule with poor (no peaks matched a fragment)
                // or missing spectra.
            }

            if (itdata->getGroup() == validation_group)
                continue;

            MolData *moldata = &(*itdata);

            // Compute the transition probabilities
            computeThetas(moldata);
            moldata->computeTransitionProbabilities();

            // Apply the peak evidence, compute the beliefs and record the sufficient
            // statistics
            if (cfg->use_single_energy_cfm) {
                beliefs_t beliefs;
                Inference infer(moldata, cfg);
                infer.calculateBeliefs(beliefs);
                recordSufficientStatistics(suft, molidx, moldata, &beliefs);
            } else {
                IPFP ipfp(moldata, cfg);
                beliefs_t *beliefs = ipfp.calculateBeliefs();
                int status = ipfp.status;
                if (status == NON_CONVERGE || status == OSC_CONVERGE)
                    num_nonconverged++;
                else if (status == COMPLETE_CONVERGE || status == CONVERGE_AFTER_MOD)
                    num_converged++;
                recordSufficientStatistics(suft, molidx, moldata, beliefs);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait
        after = time(nullptr);
        if (!cfg->use_single_energy_cfm) {
            total_numnonc = comm->collectSumInMaster(num_nonconverged);
            tot_numc = comm->collectSumInMaster(num_converged);
            std::string cvg_msg =
                    "Num Converged: " + boost::lexical_cast<std::string>(tot_numc);
            std::string noncvg_msg = "Num Non-Converged: " +
                                     boost::lexical_cast<std::string>(total_numnonc);
            if (comm->isMaster()) {
                writeStatus(cvg_msg.c_str());
                writeStatus(noncvg_msg.c_str());
            }
            comm->printToMasterOnly(cvg_msg.c_str());
            comm->printToMasterOnly(noncvg_msg.c_str());
        }
        std::string estep_time_msg =
                "[E-Step]Completed E-step processing: Time Elapsed = " +
                boost::lexical_cast<std::string>(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(estep_time_msg.c_str());
        comm->printToMasterOnly(estep_time_msg.c_str());

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master
        // Find a new set of parameters to maximize the expected log likelihood
        // (M-step)

        before = time(nullptr);
        if (cfg->ga_method == USE_LBFGS_FOR_GA)
            Q = updateParametersLBFGS(data, suft);
        else
            Q = updateParametersGradientAscent(data, suft);

        after = time(nullptr);
        std::string param_update_time_msg =
                "[M-Step]Completed M-step param update: Time Elapsed = " +
                boost::lexical_cast<std::string>(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(param_update_time_msg.c_str());
        comm->printToMasterOnly(param_update_time_msg.c_str());

        // Write the params
        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master

        before = time(nullptr);
        // validation Q value
        double valQ = 0.0;
        // Compute the final Q (with all molecules, in case only some were used in
        // the mini-batch)
        if (cfg->ga_minibatch_nth_size > 1) {
            Q = 0;
            if (comm->isMaster())
                std::cout << "Using minibatch, compute the final Q" << std::endl;
        }
        int molidx = 0, numvalmols = 0, numnonvalmols = 0;
        for (itdata = data.begin(); itdata != data.end(); ++itdata, molidx++) {
            if (itdata->getGroup() == validation_group) {
                //valQ += computeQ( molidx, *itdata, suft );
                numvalmols++;
            } else if (cfg->ga_minibatch_nth_size > 1) {
                Q += computeQ(molidx, *itdata, suft);
                numnonvalmols++;
            } else
                numnonvalmols++;
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads wait for master
        after = time(nullptr);
        std::string q_time_msg =
                "[M-Step] Finished Q compute: Time Elapsed = " +
                boost::lexical_cast<std::string>(after - before) + " seconds";
        if (comm->isMaster())
            writeStatus(q_time_msg.c_str());

        valQ = comm->collectQInMaster(valQ);
        if (cfg->ga_minibatch_nth_size > 1) {
            if (comm->isMaster()) {
                Q += addRegularizersAndUpdateGradient(nullptr);
            }
            Q = comm->collectQInMaster(Q);
            Q = comm->broadcastQ(Q);
        }
        numvalmols = comm->collectSumInMaster(numvalmols);
        numnonvalmols = comm->collectSumInMaster(numnonvalmols);

        // Check for convergence
        double Qratio = fabs((Q - prevQ) / Q);
        double bestQRatio = fabs((Q - bestQ) / Q);
        if (comm->isMaster()) {
            std::string qdif_str = "[M-Step] Q_ratio= " + boost::lexical_cast<std::string>(Qratio) + " prev_Q=" +
                                   boost::lexical_cast<std::string>(prevQ) + "\n";
            qdif_str = "Best_Q_ratio= " + boost::lexical_cast<std::string>(bestQRatio) + " best_Q=" +
                       boost::lexical_cast<std::string>(bestQ) + "\n";

            qdif_str += "Q=" + boost::lexical_cast<std::string>(Q) +
                        " ValidationQ=" + boost::lexical_cast<std::string>(valQ) + " ";
            qdif_str +=
                    "Q_avg=" + boost::lexical_cast<std::string>(Q / numnonvalmols) +
                    " ValidationQ_avg=" + boost::lexical_cast<std::string>(valQ / numvalmols) + " ";
            writeStatus(qdif_str.c_str());
            comm->printToMasterOnly(qdif_str.c_str());
        }

        // check if make any progress yet
        // two conditions: 1. Qratio is less than 1e-15
        //                 2, Q has not improved compare to the best value so far
        const double ratio_cutoff = 1e-15;
        if (bestQRatio < ratio_cutoff || bestQ > Q) {
            count_no_progress += 1;
        } // write param to file if current Q is better
        else {
            bestQ = Q;
            // Write the params
            if (comm->isMaster()) {
                std::string progress_str = "[M-Step] Found Better Q: "
                                           + boost::lexical_cast<std::string>(bestQ) + " Write to File";
                comm->printToMasterOnly(progress_str.c_str());
                writeParamsToFile(iter_out_param_filename);
                writeParamsToFile(out_param_filename);
            }
            count_no_progress = 0;
        }

        prevQ = Q;
        // check if EM meet halt flag
        if (Qratio < cfg->em_converge_thresh || count_no_progress >= 3) {
            comm->printToMasterOnly(("EM Converged after " +
                                     boost::lexical_cast<std::string>(iter) +
                                     " iterations")
                                            .c_str());
            break;
        }
        iter++;
    }

    if (iter >= MAX_EM_ITERATIONS)
        comm->printToMasterOnly(("Warning: EM did not converge after " +
                                 boost::lexical_cast<std::string>(iter) +
                                 " iterations.")
                                        .c_str());

    return bestQ;
}

void EM::initSuft(suft_counts_t &suft, std::vector<MolData> &data) {

    // Resize the suft structure for each molecule
    unsigned int num_mols = data.size();
    suft.values.resize(num_mols);
    for (unsigned int i = 0; i < num_mols; i++) {
        const FragmentGraph *fg = data[i].getFragmentGraph();
        int len = fg->getNumTransitions() + fg->getNumFragments();
        int num_spectra = data[i].getNumSpectra();
        suft.values[i].resize(len * num_spectra);
    }
}

void EM::recordSufficientStatistics(suft_counts_t &suft, int molidx,
                                    MolData *moldata, beliefs_t *beliefs) {

    const FragmentGraph *fg = moldata->getFragmentGraph();

    unsigned int num_transitions = fg->getNumTransitions();
    unsigned int num_fragments = fg->getNumFragments();

    int len_offset = num_transitions + num_fragments;

    // Accumulate the Sufficient Statistics
    for (unsigned int i = 0; i < num_transitions; i++) {

        const Transition *t = fg->getTransitionAtIdx(i);

        double belief = 0.0;
        int energy = cfg->map_d_to_energy[0];
        if (t->getFromId() == 0) // main ion is always id = 0
            belief += exp(beliefs->tn[i][0]);

        for (unsigned int d = 1; d < cfg->model_depth; d++) {
            energy = cfg->map_d_to_energy[d];
            if (energy != cfg->map_d_to_energy[d - 1]) {
                suft.values[molidx][i + cfg->map_d_to_energy[d - 1] * len_offset] =
                        belief;
                belief = 0.0;
            }
            belief += exp(beliefs->tn[i][d]);
        }
        suft.values[molidx][i + energy * len_offset] = belief;
    }

    // Accumulate the persistence terms
    int offset = num_transitions;
    for (unsigned int i = 0; i < num_fragments; i++) {

        double belief = 0.0;
        int energy = cfg->map_d_to_energy[0];

        if (i == 0) // main ion is always id = 0
            belief += exp(beliefs->ps[i][0]);
        for (unsigned int d = 1; d < cfg->model_depth; d++) {
            energy = cfg->map_d_to_energy[d];
            if (energy != cfg->map_d_to_energy[d - 1]) {
                suft.values[molidx][i + offset +
                                    cfg->map_d_to_energy[d - 1] * len_offset] = belief;
                belief = 0.0;
            }
            belief += exp(beliefs->ps[i][d]);
        }
        // TODO FIND A BETTER WAY THIS IS A SUPER HACKY FIX
        if (boost::math::isinf(belief))
            belief = 1000000000;
        suft.values[molidx][i + offset + energy * len_offset] = belief;
    }
}

static lbfgsfloatval_t lbfgs_evaluate(void *instance, const lbfgsfloatval_t *x,
                                      lbfgsfloatval_t *g, const int n,
                                      const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = ((EM *) instance)->evaluateLBFGS(x, g, n, step);
    return fx;
}

static int lbfgs_progress(void *instance, const lbfgsfloatval_t *x,
                          const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
                          const lbfgsfloatval_t xnorm,
                          const lbfgsfloatval_t gnorm,
                          const lbfgsfloatval_t step, int n, int k, int ls) {
    ((EM *) instance)->progressLBFGS(x, g, fx, xnorm, gnorm, step, n, k, ls);
    return 0;
}

double EM::updateParametersLBFGS(std::vector<MolData> &data,
                                 suft_counts_t &suft) {

    double Q = 0.0;
    lbfgs_parameter_t lparam;
    lbfgs_parameter_init(&lparam);
    lparam.delta = cfg->ga_converge_thresh;
    lparam.past = 1;
    lparam.max_iterations = cfg->ga_max_iterations;

    // Initial Q and gradient calculation (to determine used indexes - if we
    // already have them don't bother)
    if (comm->used_idxs.empty()) {
        std::vector<double> grads(param->getNumWeights(), 0.0);
        auto itdata = data.begin();
        for (int molidx = 0; itdata != data.end(); ++itdata, molidx++) {
            if (itdata->getGroup() != validation_group) {
                Q += computeAndAccumulateGradient(&grads[0], molidx, *itdata, suft,
                                                  true, comm->used_idxs);
            }
        }
        // Collect the used_idxs from all the processors into the MASTER
        comm->setMasterUsedIdxs();

        // Copy the used parameters into the LBFGS array
        if (comm->isMaster())
            zeroUnusedParams();
    }

    int N = 0;
    if (comm->isMaster())
        N = ((MasterComms *) comm)->master_used_idxs.size();
    N = comm->broadcastNumUsed(N);
    if (N > 0.1 * param->getNumWeights())
        sparse_params = false;
    lbfgsfloatval_t *x = convertCurrentParamsToLBFGS(N);
    if (!sparse_params)
        N = param->getNumWeights();

    // Select molecules to include in gradient mini-batch (will select a new set
    // at each LBFGS iteration, but for every evaluation).
    tmp_minibatch_flags.resize(data.size());
    std::vector<MolData>::iterator itdata = data.begin();
    for (int molidx = 0; itdata != data.end(); ++itdata, molidx++)
        tmp_minibatch_flags[molidx] =
                (itdata->getGroup() !=
                 validation_group); // Don't include validation molecules
    if (cfg->ga_minibatch_nth_size > 1) {
        comm->printToMasterOnly("Selecting MiniBatch for LBFGS");
        selectMiniBatch(tmp_minibatch_flags);
    }

    // Run LBFGS
    tmp_moldata_ptr_lbfgs = &data;
    tmp_suft_ptr_lbfgs = &suft;
    lbfgsfloatval_t fx;
    comm->printToMasterOnly("Running LBFGS...");
    lbfgs(N, x, &fx, lbfgs_evaluate, lbfgs_progress, this, &lparam);

    // Master converts and broadcasts final param weights and Q to all
    copyLBFGSToParams(x);
    if (sparse_params)
        lbfgs_free(x);
    Q = -fx;
    comm->broadcastParams(param.get());
    Q = comm->broadcastQ(Q);

    return (Q);
}

lbfgsfloatval_t *EM::convertCurrentParamsToLBFGS(int N) {

    lbfgsfloatval_t *x;
    if (sparse_params) {
        x = lbfgs_malloc(N);
        if (x == nullptr) {
            std::cout << "ERROR: Failed to allocate a memory block for variables."
                      << std::endl;
            throw EMComputationException();
        }

        // Only fill the actual parameters in the master (the others we can update
        // at the start of  each evaluate call).
        if (comm->isMaster()) {
            std::set<unsigned int>::iterator it =
                    ((MasterComms *) comm)->master_used_idxs.begin();
            for (unsigned int i = 0;
                 it != ((MasterComms *) comm)->master_used_idxs.end(); ++it)
                x[i++] = param->getWeightAtIdx(*it);
        }
    } else
        x = &((*param->getWeightsPtr())[0]); // If not sparse, just use the existing
    // param array
    return x;
}

void EM::copyGradsToLBFGS(lbfgsfloatval_t *g, std::vector<double> &grads,
                          int n) {

    if (comm->isMaster()) {
        if (sparse_params) {
            std::set<unsigned int>::iterator it =
                    ((MasterComms *) comm)->master_used_idxs.begin();
            for (unsigned int i = 0;
                 it != ((MasterComms *) comm)->master_used_idxs.end(); ++it)
                g[i++] = -grads[*it];
        } else {
            for (int i = 0; i < n; i++)
                g[i] *= -1;
        }
    }
    comm->broadcastGorX(g, n);
}

void EM::copyLBFGSToParams(const lbfgsfloatval_t *x) {
    if (sparse_params && comm->isMaster()) {
        std::set<unsigned int>::iterator it =
                ((MasterComms *) comm)->master_used_idxs.begin();
        for (unsigned int i = 0;
             it != ((MasterComms *) comm)->master_used_idxs.end(); ++it)
            param->setWeightAtIdx(x[i++], *it);
    }
}

lbfgsfloatval_t EM::evaluateLBFGS(const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
                                  const int n, const lbfgsfloatval_t step) {

    // Retrieve params from LBFGS and sync between processors
    copyLBFGSToParams(x);
    comm->broadcastParams(param.get());

    // Set the location for the computed gradients (new array if sparse, else the
    // provided g) and initialise
    double *grad_ptr;
    std::vector<double> grads;
    if (sparse_params) {
        grads.resize(param->getNumWeights(), 0.0);
        grad_ptr = &grads[0];
    } else {
        grad_ptr = g;
        std::set<unsigned int>::iterator sit = comm->used_idxs.begin();
        for (; sit != comm->used_idxs.end(); ++sit)
            *(grad_ptr + *sit) = 0.0;
    }

    // Compute Q and the gradient
    double Q = 0.0;
    auto itdata = tmp_moldata_ptr_lbfgs->begin();
    for (int molidx = 0; itdata != tmp_moldata_ptr_lbfgs->end();
         ++itdata, molidx++) {
        if (tmp_minibatch_flags[molidx])
            Q += computeAndAccumulateGradient(grad_ptr, molidx, *itdata,
                                              *tmp_suft_ptr_lbfgs, false,
                                              comm->used_idxs);
    }

    if (comm->isMaster())
        Q += addRegularizersAndUpdateGradient(grad_ptr);
    comm->collectGradsInMaster(grad_ptr);
    Q = comm->collectQInMaster(Q);
    Q = comm->broadcastQ(Q);

    // Move the computed gradients into the lbfgs structure (note: only used idxs
    // are included)
    copyGradsToLBFGS(g, grads, n);
    return -Q;
}

void EM::progressLBFGS(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g,
                       const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm,
                       const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step,
                       int n, int k, int ls) {

    if (comm->isMaster()) {
        writeStatus(("LBFGS Iteration " + boost::lexical_cast<std::string>(k) +
                     ": fx = " + boost::lexical_cast<std::string>(fx))
                            .c_str());
        std::cout << "LBFGS Iteration " << k << ": fx = " << fx << std::endl;
    }

    // Select molecules to include in next gradient mini-batch (will select a new
    // set at each LBFGS iteration, but not for every evaluation).
    std::vector<MolData>::iterator itdata = tmp_moldata_ptr_lbfgs->begin();
    for (int molidx = 0; itdata != tmp_moldata_ptr_lbfgs->end();
         ++itdata, molidx++)
        tmp_minibatch_flags[molidx] =
                (itdata->getGroup() !=
                 validation_group); // Don't include validation molecules
    if (cfg->ga_minibatch_nth_size > 1)
        selectMiniBatch(tmp_minibatch_flags);
}

double EM::updateParametersGradientAscent(std::vector<MolData> &data,
                                          suft_counts_t &suft) {

    // DBL_MIN is the smallest positive double
    // -DBL_MAX is the smallest negative double
    double Q = 0.0, prevQ = -DBL_MAX, bestQ = -DBL_MAX;

    std::vector<double> grads(param->getNumWeights(), 0.0);
    std::vector<double> prev_v(param->getNumWeights(), 0.0);

    // for adam and AMSgrad
    std::vector<double> first_moment_vector(param->getNumWeights(), 0.0);
    std::vector<double> second_moment_vector(param->getNumWeights(), 0.0);
    std::vector<double> second_moment_max_vector(param->getNumWeights(), 0.0);

    // for adaDelta
    std::vector<double> mean_squared_gradients(param->getNumWeights(), 0.0);
    std::vector<double> mean_squared_delta_x(param->getNumWeights(), 0.0);

    // Initial Q and gradient calculation (to determine used indexes)
    if (comm->used_idxs.empty()) {
        auto itdata = data.begin();
        for (int molidx = 0; itdata != data.end(); ++itdata, molidx++) {
            if (itdata->getGroup() != validation_group) {
                Q += computeAndAccumulateGradient(&grads[0], molidx, *itdata, suft,
                                                  true, comm->used_idxs);
            }
        }
        comm->setMasterUsedIdxs();
        if (comm->isMaster())
            zeroUnusedParams();
    }
    int N = 0;
    if (comm->isMaster())
        N = ((MasterComms *) comm)->master_used_idxs.size();
    N = comm->broadcastNumUsed(N);

    int iter = 0;
    double learn_mult = 1.0;

    //TODO maybe I should force iteration time such time ga covers all data points
    int max_iteration = cfg->ga_max_iterations;
    int no_progress_count = 0;
    while (iter++ < max_iteration
           && fabs((Q - prevQ) / Q) >= cfg->ga_converge_thresh
           && no_progress_count < 3) {

        if (Q < prevQ && iter > 1 && cfg->ga_method == USE_MOMENTUM_FOR_GA)
            learn_mult = learn_mult * 0.5;

        // adjust learning rate
        double learn_rate = cfg->starting_step_size * learn_mult;
        if (USE_DEFAULT_DECAY == cfg->ga_decay_method)
            learn_rate *= 1.0 / (1.0 + cfg->decay_rate * (iter - 1));
        else if (USE_EXP_DECAY == cfg->ga_decay_method)
            learn_rate *= std::exp(-cfg->exp_decay_k * iter);
        else if (USE_STEP_DECAY == cfg->ga_decay_method)
            learn_rate *= std::pow(cfg->step_decay_drop, std::floor(iter / cfg->step_decay_epochs_drop));

        if (iter > 1)
            prevQ = Q;

        // Select molecules to include in gradient mini-batch.
        std::vector<int> minibatch_flags(data.size());
        auto itdata = data.begin();
        for (int molidx = 0; itdata != data.end(); ++itdata, molidx++) {
            // Don't include validation molecules
            minibatch_flags[molidx] = (itdata->getGroup() != validation_group);
        }

        if (cfg->ga_minibatch_nth_size > 1) {
            selectMiniBatch(minibatch_flags);
        }

        // Compute Q and the gradient
        std::fill(grads.begin(), grads.end(), 0.0);
        Q = 0.0;
        itdata = data.begin();
        for (int molidx = 0; itdata != data.end(); ++itdata, molidx++) {
            if (minibatch_flags[molidx])
                Q += computeAndAccumulateGradient(&grads[0], molidx, *itdata, suft,
                                                  false, comm->used_idxs);
        }
        if (comm->isMaster())
            Q += addRegularizersAndUpdateGradient(&grads[0]);
        comm->collectGradsInMaster(&grads[0]);

        Q = comm->collectQInMaster(Q);
        Q = comm->broadcastQ(Q);

        if (comm->isMaster())
            std::cout << iter << ":  Q=" << Q << " prevQ=" << prevQ << " Learning_Rate= " << learn_rate
                      << std::endl;

        // Step the parameters
        if (comm->isMaster()) {
            switch (cfg->ga_method) {
                case USE_MOMENTUM_FOR_GA:
                    param->adjustWeightsByGrads_Momentum(grads,
                                                         ((MasterComms *) comm)->master_used_idxs,
                                                         learn_rate, cfg->ga_momentum, prev_v);
                    break;
                case USE_ADAM_FOR_GA:
                    param->adjustWeightsByGrads_Adam(grads,
                                                     ((MasterComms *) comm)->master_used_idxs,
                                                     learn_rate,
                                                     cfg->ga_adam_beta_1,
                                                     cfg->ga_adam_beta_2,
                                                     cfg->ga_eps,
                                                     iter,
                                                     first_moment_vector,
                                                     second_moment_vector);
                    break;
                case USE_AMSGRAD_FOR_GA:
                    param->adjustWeightsByGrads_AMSgrad(grads,
                                                        ((MasterComms *) comm)->master_used_idxs,
                                                        learn_rate,
                                                        cfg->ga_adam_beta_1,
                                                        cfg->ga_adam_beta_2,
                                                        cfg->ga_eps,
                                                        iter,
                                                        first_moment_vector,
                                                        second_moment_vector,
                                                        second_moment_max_vector);
                    break;
                case USE_ADADELTA_FOR_GA:
                    param->adjustWeightsByGrads_Adadelta(grads,
                                                         ((MasterComms *) comm)->master_used_idxs,
                                                         learn_rate,
                                                         cfg->ga_adadelta_rho,
                                                         cfg->ga_eps,
                                                         mean_squared_gradients,
                                                         mean_squared_delta_x);
                    break;
                default:
                    param->adjustWeightsByGrads_Momentum(grads,
                                                         ((MasterComms *) comm)->master_used_idxs,
                                                         learn_rate, cfg->ga_momentum, prev_v);
            }
        }
        comm->broadcastParams(param.get());

        if (cfg->ga_use_best_q) {
            if (bestQ > Q) {
                no_progress_count++;
            } else {
                bestQ = Q;
                no_progress_count = 0;
            }
        } else {
            if (prevQ > Q) {
                no_progress_count++;
            } else {
                no_progress_count = 0;
            }
        }

    }

    if (comm->isMaster()) {
        if (iter == cfg->ga_max_iterations)
            std::cout << "Gradient ascent did not converge" << std::endl;
        else
            std::cout << "Gradient ascent converged after " << iter << " iterations"
                      << std::endl;
    }
    return Q;
}

double EM::computeAndAccumulateGradient(double *grads, int molidx,
                                        MolData &moldata, suft_counts_t &suft,
                                        bool record_used_idxs_only,
                                        std::set<unsigned int> &used_idxs) {

    double Q = 0.0;
    const FragmentGraph *fg = moldata.getFragmentGraph();
    unsigned int num_transitions = fg->getNumTransitions();
    unsigned int num_fragments = fg->getNumFragments();

    int offset = num_transitions;

    if (!moldata.hasComputedGraph())
        return Q;

    // Compute the latest transition thetas
    if (!record_used_idxs_only)
        moldata.computeTransitionThetas(*param);
    suft_t *suft_values = &(suft.values[molidx]);

    // Collect energies to compute
    std::vector<unsigned int> energies;
    getEnergiesLevels(energies);

    // Compute the gradients
    for (auto energy : energies) {

        unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();
        unsigned int suft_offset = energy * (num_transitions + num_fragments);

        std::vector<int> not_so_random_selected;
        // NOTE Make sure computeTransitionThetas before this
        // otherwise this will crash
        if (USE_GRAPH_RANDOM_SAMPLING == cfg->ga_sampling_method && !record_used_idxs_only) {
            moldata.getSampledTransitionIds(not_so_random_selected,
                                            cfg->ga_graph_sampling_k,
                                            energy, m_rng,
                                            m_uniform_dist);
            std::sort(not_so_random_selected.begin(), not_so_random_selected.end());

        }

        // Iterate over from_id (i)
        auto frag_trans_map = fg->getFromIdTMap()->begin();

        for (int from_idx = 0; frag_trans_map != fg->getFromIdTMap()->end(); ++frag_trans_map, from_idx++) {

            if (record_used_idxs_only) {
                for (auto trans_id : *frag_trans_map) {
                    const FeatureVector *fv = moldata.getFeatureVectorForIdx(trans_id);
                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                        used_idxs.insert(fv_it->first + grad_offset);
                    }
                }
            } else {
                //Do some random selection
                std::vector<int> frag_trans_local_copy;
                for (auto item: *frag_trans_map) {
                    frag_trans_local_copy.emplace_back(item);
                }
                if (USE_RANDOM_SAMPLING == cfg->ga_sampling_method) {
                    std::shuffle(frag_trans_local_copy.begin(), frag_trans_local_copy.end(), m_rng);
                    auto resize_len = (unsigned) ((double) frag_trans_local_copy.size() *
                                                  cfg->random_sampling_threshold);
                    if (resize_len == 0)
                        resize_len = 1;
                    frag_trans_local_copy.resize(resize_len);
                } else if (USE_GRAPH_RANDOM_SAMPLING == cfg->ga_sampling_method) {
                    if (!frag_trans_local_copy.empty()) {
                        std::vector<int> v_intersection;
                        std::sort(frag_trans_local_copy.begin(), frag_trans_local_copy.end());
                        std::set_intersection(frag_trans_local_copy.begin(), frag_trans_local_copy.end(),
                                              not_so_random_selected.begin(), not_so_random_selected.end(),
                                              std::back_inserter(v_intersection));
                        frag_trans_local_copy = v_intersection;
                    }
                }


                // Calculate the denominator of the sum terms
                double denom = 1.0;
                for (auto trans_id : frag_trans_local_copy)
                    denom += exp(moldata.getThetaForIdx(energy, trans_id));

                // Complete the innermost sum terms	(sum over j')
                std::map<unsigned int, double> sum_terms;

                for (auto trans_id : frag_trans_local_copy) {
                    const FeatureVector *fv = moldata.getFeatureVectorForIdx(trans_id);

                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                        auto fv_idx = fv_it->first;
                        double val = exp(moldata.getThetaForIdx(energy, trans_id)) / denom;
                        if (sum_terms.find(fv_idx) != sum_terms.end())
                            sum_terms[fv_idx] += val;
                        else
                            sum_terms[fv_idx] = val;
                    }

                }

                // Accumulate the transition (i \neq j) terms of the gradient (sum over j)
                double nu_sum = 0.0;
                for (auto trans_id : frag_trans_local_copy) {
                    double nu = (*suft_values)[trans_id + suft_offset];
                    nu_sum += nu;
                    const FeatureVector *fv = moldata.getFeatureVectorForIdx(trans_id);

                    for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
                        auto fv_idx = fv_it->first;
                        *(grads + fv_idx + grad_offset) += nu;
                    }

                    Q += nu * (moldata.getThetaForIdx(energy, trans_id) - log(denom));
                    if (boost::math::isnan(Q))
                        std::cerr << "Setp1" << std::endl;
                }



                // Accumulate the last term of each transition and the
                // persistence (i = j) terms of the gradient and Q
                double nu = (*suft_values)[offset + from_idx + suft_offset]; // persistence (i=j)
                for (auto &sit: sum_terms) {
                    *(grads + sit.first + grad_offset) -= (nu_sum + nu) * sit.second;
                }
                Q -= nu * log(denom);
                if (boost::math::isnan(Q)) {
                    std::cerr << moldata.getId() << "Setp 2" << nu << " " << denom << " " << from_idx + suft_offset
                              << std::endl;
                    break;
                }
            }
        }
    }
    if (boost::math::isnan(Q))
        std::cerr << moldata.getId() << std::endl;
    return Q;
}

double EM::computeQ(int molidx, MolData &moldata, suft_counts_t &suft) {

    double Q = 0.0;
    const FragmentGraph *fg = moldata.getFragmentGraph();
    unsigned int num_transitions = fg->getNumTransitions();
    unsigned int num_fragments = fg->getNumFragments();

    int offset = num_transitions;

    if (!moldata.hasComputedGraph())
        return Q;

    // Compute the latest transition thetas
    moldata.computeTransitionThetas(*param);
    suft_t *suft_values = &(suft.values[molidx]);

    // Collect energies to compute

    std::vector<unsigned int> energies;
    getEnergiesLevels(energies);

    // Compute
    for (auto energy : energies) {
        unsigned int suft_offset = energy * (num_transitions + num_fragments);

        // Iterate over from_id (i)
        auto it = fg->getFromIdTMap()->begin();
        for (int from_idx = 0; it != fg->getFromIdTMap()->end(); ++it, from_idx++) {

            // Calculate the denominator of the sum terms
            double denom = 1.0;
            for (auto itt : *it)
                denom += exp(moldata.getThetaForIdx(energy, itt));

            // Accumulate the transition (i \neq j) terms of the gradient (sum over j)
            double nu_sum = 0.0;
            for (auto itt : *it) {
                double nu = (*suft_values)[itt + suft_offset];
                Q += nu * (moldata.getThetaForIdx(energy, itt) - log(denom));
            }

            // Accumulate the last term of each transition and the
            // persistence (i = j) terms of the gradient and Q
            double nu =
                    (*suft_values)[offset + from_idx + suft_offset]; // persistence (i=j)
            Q -= nu * log(denom);
        }
    }
    return Q;
}

double EM::addRegularizersAndUpdateGradient(double *grads) {

    double Q = 0.0;
    auto it = ((MasterComms *) comm)->master_used_idxs.begin();
    for (; it != ((MasterComms *) comm)->master_used_idxs.end(); ++it) {

        double weight = param->getWeightAtIdx(*it);
        Q -= 0.5 * cfg->lambda * weight * weight;
        if (grads != nullptr)
            *(grads + *it) -= cfg->lambda * weight;
    }

    // Remove the Bias terms (don't regularize the bias terms!)
    unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
    for (unsigned int energy = 0; energy < param->getNumEnergyLevels();
         energy++) {
        double bias = param->getWeightAtIdx(energy * weights_per_energy);
        Q += 0.5 * cfg->lambda * bias * bias;
        if (grads != nullptr)
            *(grads + energy * weights_per_energy) += cfg->lambda * bias;
    }
    return Q;
}


void EM::zeroUnusedParams() {

    unsigned int i;
    for (i = 0; i < param->getNumWeights(); i++) {
        if (((MasterComms *) comm)->master_used_idxs.find(i) ==
            ((MasterComms *) comm)->master_used_idxs.end())
            param->setWeightAtIdx(0.0, i);
    }
}

void EM::getEnergiesLevels(std::vector<unsigned int> &energies) {
    unsigned int energy;
    int prev_energy = -1;
    for (unsigned int d = 0; d < cfg->model_depth; d++) {
        energy = cfg->map_d_to_energy[d];
        if (energy != prev_energy)
            energies.push_back(energy);
        prev_energy = energy;
    }
}


void EM::selectMiniBatch(std::vector<int> &initialized_minibatch_flags) {

    // The flags are initialized to 1's for selectable molecules, so set unwanted
    // molecule flags to 0
    int num_mols = initialized_minibatch_flags.size();
    std::vector<int> idxs(num_mols);
    int count = 0;
    for (int i = 0; i < num_mols; i++)
        if (initialized_minibatch_flags[i])
            idxs[count++] = i;
    idxs.resize(count);


    std::shuffle(idxs.begin(), idxs.end(), m_rng);

    int num_minibatch_mols =
            (num_mols + cfg->ga_minibatch_nth_size - 1) / cfg->ga_minibatch_nth_size;
    for (int i = num_minibatch_mols; i < idxs.size(); i++)
        initialized_minibatch_flags[idxs[i]] = 0;
}