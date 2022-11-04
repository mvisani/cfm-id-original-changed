/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FeatureCalculator.h
#
# Description: 	Code for computing features for fragmentations.
#
# Copyright (c) 2018
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#pragma once

#include "Util.h"
#include "FunctionalGroups.h"
#include "Feature.h"
#include "FeatureVector.h"

#include <GraphMol/FragCatalog/FragCatParams.h>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/lexical_cast.hpp>

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

// Exception to throw when the input feature configuration file is invalid
class InvalidConfigException : public std::exception {

    virtual const char *what() const noexcept {
        return "Invalid Feature Configuration File";
    }
};

class FeatureCalculationException : public std::exception {
private:
    std::string message_;

public:
    FeatureCalculationException(const std::string &message) noexcept
            : message_(message) {};

    virtual const char *what() const noexcept {
        std::cout << "Error computing feature vector: " << message_ << std::endl;
        return message_.c_str();
    }

    ~FeatureCalculationException() noexcept {};
};

// Class to compute a feature vector
class FeatureCalculator {

public:
    // Constructor: Initialise the calculator using a config file listing features
    FeatureCalculator(std::string &config_filename);

    // Constructor: Initialise the calculator using a list of feature names
    FeatureCalculator(std::vector<std::string> &feature_list);

    // Compute the expected number of total features
    unsigned int getNumFeatures();

    // Retrieve the list of feature names being used
    std::vector<std::string> getFeatureNames();

    // Retrieve a list of valid feature names (for testing)
    static const std::vector<std::string> getValidFeatureNames();

    // Compute the feature vector for the input ion and nl (with labeled Root
    // atoms)
    // - NB: responsibility of caller to delete.
    FeatureVector *computeFeatureVector(const RootedROMol *ion, const RootedROMol *nl, const romol_ptr_t precursor_ion);

    bool includesFeature(const std::string &fname);

private:
    // List of feature classes ready to be used
    static const boost::ptr_vector<BreakFeature> &breakFeatureCogs();

    // List of fragmentation features ready to be used
    static const boost::ptr_vector<FragmentFeature> &fragmentFeatureCogs();

    // Indexes of feature classes that are selected for use
    std::vector<int> used_break_feature_idxs;
    std::vector<int> used_fragement_feature_idxs;

    // Helper function - Configure feature for use
    void configureFeature(std::string &name);
};

