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

#include "Param.h"
#include "NNParam.h"
#include <string>
#include <set>

static const int MASTER = 0;

class Comms {

public:
    Comms();

    virtual void setMasterUsedIdxs() = 0;

    virtual void collectGradsInMaster(double *grads) = 0;

    virtual void broadcastParams(Param *param) = 0;

    virtual void printToMasterOnly(const char *msg) = 0;

    void broadcastParamsOrigMpi(Param *param);

    void printWithWorkerId(const char *msg);

    double collectQInMaster(double Q);

    bool isMaster() { return mpi_rank == MASTER; };

    int broadcastConverged(int converged);

    int broadcastNumUsed(int num_used);

    bool broadcastBooleanFlag(bool flag);

    int collectSumInMaster(int partial);

    double broadcastQ(double Q);

    std::set<unsigned int> used_idxs;
    unsigned int num_used;

    virtual ~Comms() = default;

protected:
    int mpi_rank;
    int mpi_nump;
};

class WorkerComms : public Comms {

public:
    void setMasterUsedIdxs() override;

    void collectGradsInMaster(double *grads) override;

    void broadcastParams(Param *param) override;

    void printToMasterOnly(const char *msg) override {};    //Do nothing
};

class MasterComms : public Comms {

public:
    void setMasterUsedIdxs() override;

    void collectGradsInMaster(double *grads) override;

    void broadcastParams(Param *param) override;

    std::set<unsigned int> master_used_idxs;
    std::vector<std::set<unsigned int> > worker_used_idxs;
    std::vector<unsigned int> worker_num_used;

    void printToMasterOnly(const char *msg) override;
};

#endif // __COMMS_H__
