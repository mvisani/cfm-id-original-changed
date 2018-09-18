/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# nn_param.h
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

#ifndef __NN_PARAM_H__
#define __NN_PARAM_H__


#include "Param.h"

#include <string>

typedef std::vector<double> azd_vals_t;


//Exception to throw when the activation function id is unknown
class NNParamActivationFunctionIdException : public std::exception {

    virtual const char *what() const noexcept {
        return "Unknown activation function";
    }
};

//Exception to throw when a nn_param file is expected to contain the neural net configuration
//information but doesn't
class NNParamFileReadException : public std::exception {

    virtual const char *what() const noexcept {
        return "Couldn't find neural net configuration information in nn_param file";
    }
};


class NNParam : public Param {
public:
    NNParam(std::vector<std::string> a_feature_list, int a_num_energy_levels,
                std::vector<int> &a_hlayer_num_nodes, std::vector<int> &a_act_func_ids,
                std::vector<double> &a_dropout_probs);

    //Constructor for loading parameters from file
    NNParam(std::string &filename);

    //Save parameters to file (non-sparse format, includes neural net configuration)
    void saveToFile(std::string &filename);

    //Compute the theta value for an input feature vector and energy
    //based on the current weight settings
    //is should only be used in prediction phase or compute loss
    double computeTheta(const FeatureVector &fv, int energy);

    // Compute the theta value for an input feature vector and energy in training time
    // this function is added to use drop out in nerual network
    // which in Inverted Dropout Implementation, drop out only chance FowardPassing phase
    // this function also keeps a record of the z and a values along the way
    // (for forwards step in forwards-backwards algorithm)
    // use_dropout flag should set to false during used_idx collection
    double computeTheta(const FeatureVector &fv, int energy, azd_vals_t &z_values, azd_vals_t &a_values,
                        bool already_sized = false, bool use_dropout = false);

    void
    computeDeltas(std::vector<azd_vals_t> &deltasA, std::vector<azd_vals_t> &deltasB, std::vector<azd_vals_t> &z_values,
                  std::vector<azd_vals_t> &a_values, double rho_denom, int energy);

    void
    computeUnweightedGradients(std::vector<std::vector<double> > &unweighted_grads, std::set<unsigned int> &used_idxs,
                               std::vector<const FeatureVector *> &fvs, std::vector<azd_vals_t> &deltasA,
                               std::vector<azd_vals_t> &deltasB, std::vector<azd_vals_t> &a_values);

    unsigned int getTotalNumNodes() const { return total_nodes; };    //Only hidden and output nodes (not input nodes)
    unsigned int getSecondLayerWeightOffset() const { return hlayer_num_nodes[0] * expected_num_input_features; };

    virtual void getBiasIndexes(std::vector<unsigned int> &bias_indexes) override;

    virtual void initWeights(int init_type) override;

    virtual boost::container::vector<bool> *getDropoutsPtr() override { return &is_dropped; };

    virtual void updateDropoutsRate(double delta, double lower_bound) {
        for(auto & hlayer_dropout_prob : hlayer_dropout_probs){
            if(hlayer_dropout_prob + delta > lower_bound){
                hlayer_dropout_prob += delta;
            }
        }
    };
    
    void rollDropouts(int iter, double delta) override;

protected:
    //Initialisation options
    virtual void randomUniformInit() override;
    virtual void randomNormalInit() override;
    void varianceScalingInit();

private:
    // hold num of node for each hlayer
    std::vector<int> hlayer_num_nodes;
    // hold drop out prob for each hlayer
    std::vector<double> hlayer_dropout_probs;
    // addition paramenter for drop outs
    boost::container::vector<bool> is_dropped;

    // tmp value for  varianceScalingInit
    std::vector<int> num_weights_per_layer;

    unsigned int total_nodes;
    unsigned int input_layer_node_num;

    //Activation functions and their derivatives (kept general in case we want to try other activation functions...)
    std::vector<int> act_func_ids;
    std::vector<double (*)(double)> act_funcs;
    std::vector<double (*)(double)> deriv_funcs;

    //Linear Activation (usually used for the last layer)
    static double linear_activation(double input) { return input; };

    static double linear_derivative(double input) { return 1; };

    //ReLU Activation
    static double relu_activation(double input) { if (input > 0) return input; else return 0; };

    static double relu_derivative(double input) { if (input > 0) return 1; else return 0; };

    static double leaky_relu_activation(double input) { if (input > 0) return input; else return input * 0.01; };

    static double leaky_relu_derivative(double input) { if (input > 0) return 1; else return 0.01; };

    static double neg_relu_activation(double input) { if (input < 0) return input; else return 0; };

    static double neg_relu_derivative(double input) { if (input < 0) return 1; else return 0; };

    //Function to configure activation functions used in each layer
    void setActivationFunctionsFromIds();

    azd_vals_t tmp_z_values, tmp_a_values;

};


#endif // __NN_PARAM_H__
