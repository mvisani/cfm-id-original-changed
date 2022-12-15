/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# comms.h
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see param.cpp.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __COMMS_H__
#define __COMMS_H__

#include "mpi.h"

#include "Param.h"
#include "NNParam.h"
#include <string>
#include <set>

static const int MASTER = 0;

class Comms {

public:
    Comms();

    virtual void setMasterUsedIdxs() = 0;

    virtual void collectGradsInMaster(std::vector<float> &grads) = 0;

    void collectGradsInMasterOrigMpi(std::vector<float> &grads);

    virtual void broadcastParamsWeights(Param *param) = 0;

    virtual void printToMasterOnly(const char *msg) = 0;

    void broadcastParamsWeightsOrigMpi(Param *param);

    void broadcastDropouts(Param *param);

    void printWithWorkerId(const char *msg);

    float collectQInMaster(float Q);

    bool isMaster() { return mpi_rank == MASTER; };

    int broadcastCountValue(int value);

    int collectSumInMaster(int partial);

    float broadcastQ(float Q);

    float getTimeUsages(float time_used, MPI_Op op);

    void gatherTimeUsages(float time_used, std::vector<float> &time_used_vector);

    std::set<unsigned int> used_idxs;
    unsigned int num_used;

    int getNumProcesses() { return mpi_nump; };
    virtual ~Comms() = default;

protected:
    int mpi_rank;
    int mpi_nump;
};

class WorkerComms : public Comms {

public:
    void setMasterUsedIdxs() override;

    void collectGradsInMaster(std::vector<float> &grads) override;

    void broadcastParamsWeights(Param *param) override;

    void printToMasterOnly(const char *msg) override {};    //Do nothing
};

class MasterComms : public Comms {

public:
    void setMasterUsedIdxs() override;

    void collectGradsInMaster(std::vector<float> &grads) override;

    void broadcastParamsWeights(Param *param) override;

    std::set<unsigned int> master_used_idxs;
    std::vector<std::set<unsigned int> > worker_used_idxs;
    std::vector<unsigned int> worker_num_used;

    void printToMasterOnly(const char *msg) override;
};

#endif // __COMMS_H__
