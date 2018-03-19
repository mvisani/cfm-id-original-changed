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

//Exception to throw when the input feature vector configuration doesn't match the parameters
class ParamFeatureMismatchException : public std::exception {

    virtual const char *what() const throw() {
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
    Param(std::string &filename);

    //Append a set of parameters to the current parameters
    //Either all energy levels to the next slot (if energy < 0),
    //or just the energy level specified.
    void appendNextEnergyParams(Param &next_param, int energy = -1);

    //Copy the highest energy parameters to create another higher energy level
    void appendRepeatedPrevEnergyParams();

    //Initialisation options
    virtual void randomInit();

    void zeroInit();

    void fullZeroInit();

    //Compute the theta value for an input feature vector and energy based
    //on the current weight settings
    virtual double computeTheta(const FeatureVector &fv, int energy);

    // Update Weight use Momentum + Nesterov
    void adjustWeightsByGrads_Momentum(std::vector<double> &grads,
                                       std::set<unsigned int> &used_idxs,
                                       double learning_rate,
                                       double momentum, std::vector<double> &prev_v);
    // Update Weight use ADAM
    void adjustWeightsByGrads_Adam(std::vector<double> &grads,
                                   std::set<unsigned int> &used_idxs,
                                   const double & learning_rate,
                                   const double & beta1,
                                   const double & beta2,
                                   const double & eps,
                                   const int & iteration_count,
                                   std::vector<double> &first_moment_vector,
                                   std::vector<double> &second_moment_vector);

    // Update Weight use AdaDelta
    void adjustWeightsByGrads_Adadelta(std::vector<double> &grads,
                                       std::set<unsigned int> &used_idxs,
                                       const double & learning_rate,
                                       const double & decay_rate,
                                       const double & eps,
                                       std::vector<double> &mean_squared_gradients,
                                       std::vector<double> &mean_squared_delta_x);

    //Set the value of a weight
    void setWeightAtIdx(double value, int index) { weights[index] = value; };

    //Print the current parameters to the given output stream
    void reportParameters(std::ostream &out);

    //Save parameters to file
    virtual void saveToFile(std::string &filename);

    //Access functions
    double getWeightAtIdx(int index) { return weights[index]; };

    std::vector<double> *getWeightsPtr() { return &weights; };

    unsigned int getNumWeights() { return weights.size(); };

    unsigned int getNumWeightsPerEnergyLevel() { return weights.size() / num_energy_levels; };

    unsigned int getNumEnergyLevels() { return num_energy_levels; };

    std::vector<std::string> *getFeatureNames() { return &feature_list; };

protected:
    std::vector<double> weights;
    unsigned int num_energy_levels;
    std::vector<std::string> feature_list;
    int expected_num_input_features;

};

#endif // __PARAM_H__
