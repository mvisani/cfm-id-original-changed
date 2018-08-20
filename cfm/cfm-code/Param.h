/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# param.h
#
# Description: 	Class for parameterization of fragmentation probabilities
#				within the bayesian network fragmentation trees.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __PARAM_H__
#define __PARAM_H__

#include "Feature.h"

#include <string>
#include <boost/container/vector.hpp>

//Exception to throw when the input feature vector configuration doesn't match the parameters
class ParamFeatureMismatchException : public std::exception {

    const char *what() const noexcept override {
        return "Mismatch between feature vector length and num parameters";
    }
};


class Param {
public:
    //Constructor to initialise parameter weight size from a feature list
    Param(std::vector<std::string> a_feature_list, int a_num_energy_levels);

    //Constructor to create skeleton parameter copy
    //(doesn't actually copy weights, just resizes to same)
    Param(Param &param_instance);

    //Constructor for loading parameters from file
    explicit Param(std::string &filename);

    //Append a set of parameters to the current parameters
    //Either all energy levels to the next slot (if energy < 0),
    //or just the energy level specified.
    void appendNextEnergyParams(Param &next_param, int energy = -1);

    //Copy the highest energy parameters to create another higher energy level
    void appendRepeatedPrevEnergyParams();

    //Compute the theta value for an input feature vector and energy based
    //on the current weight settings
    virtual double computeTheta(const FeatureVector &fv, int energy);

    //Set the value of a weight
    void setWeightAtIdx(double value, int index) { weights[index] = value; };

    //Save parameters to file
    virtual void saveToFile(std::string &filename);

    //Access functions
    double getWeightAtIdx(int index) { return weights[index]; };

    std::vector<double> *getWeightsPtr() { return &weights; };

    //this will be changed once we add dropouts for linear model, for now
    //it only return nullptr
    // NOTE use boost vector of bool for mpi
    // because std vector can not return bool&
    virtual boost::container::vector<bool> *getDropoutsPtr() { return nullptr; };

    unsigned int getNumWeights() { return weights.size(); };

    unsigned int getNumWeightsPerEnergyLevel() { return weights.size() / num_energy_levels; };

    unsigned int getNumEnergyLevels() { return num_energy_levels; };

    std::vector<std::string> *getFeatureNames() { return &feature_list; };

    virtual void getBiasIndexes(std::vector<unsigned int> &bias_indexes){ bias_indexes.push_back(0);};

    virtual void initWeights(int init_type);

    // do nothing roolDropots function
    // we may add Dropouts for linear model later
    virtual void rollDropouts() {};

protected:
    std::vector<double> weights;
    unsigned int num_energy_levels;
    std::vector<std::string> feature_list;
    int expected_num_input_features;

    //Initialisation options
    virtual void randomUniformInit();
    virtual void randomNormalInit();
    void zeroInit();
    void fullZeroInit();

};

#endif // __PARAM_H__
