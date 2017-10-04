/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Features.h
#
# Description: 	Code for computing features for fragmentations.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FEATURE_H__
#define __FEATURE_H__

#include "Util.h"
#include "FunctionalGroups.h"

#include <GraphMol/FragCatalog/FragCatParams.h>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/lexical_cast.hpp>

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

struct input_file_t;

typedef std::pair<std::string, std::string> symbol_pair_t;

// Exception to throw when the input feature configuration file is invalid
class InvalidConfigException : public std::exception {

  virtual const char *what() const throw() {
    return "Invalid Feature Configuration File";
  }
};

class FeatureCalculationException : public std::exception {
private:
  std::string message_;

public:
  FeatureCalculationException(const std::string &message) throw()
      : message_(message){};
  virtual const char *what() const throw() {
    std::cout << "Error computing feature vector: " << message_ << std::endl;
    return message_.c_str();
  }
  ~FeatureCalculationException() throw(){};
};

// Structure to hold a sparse computed feature vector
typedef unsigned int feature_t;
class FeatureVector {
public:
  FeatureVector() { fv_idx = 0; };
  void addFeature(double value);
  void addFeatureAtIdx(double value, unsigned int idx);
  unsigned int getTotalLength() const { return fv_idx; };
  feature_t const getFeature(int idx) const { return fv[idx]; };
  std::vector<feature_t>::const_iterator getFeatureBegin() const {
    return fv.begin();
  };
  std::vector<feature_t>::const_iterator getFeatureEnd() const {
    return fv.end();
  };
  unsigned int getNumSetFeatures() const { return fv.size(); };

private:
  std::vector<feature_t> fv;
  unsigned int fv_idx;
};

// Base class to compute a feature - all features should inherit from this
class Feature {

public:
  virtual void compute(FeatureVector &fv, const RootedROMolPtr *ion,
                       const RootedROMolPtr *nl) const = 0;
  unsigned int getSize() const { return size; };
  std::string getName() const { return name; };
  virtual ~Feature(){};

protected:
  unsigned int size;
  std::string name;
  static const std::vector<std::string> &OKsymbols();
  static const std::vector<std::string> &OKSymbolsLess();
  void replaceUncommonWithX(std::string &symbol) const;
  int getSymbolsLessIndex(std::string &symbol) const;
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
  FeatureVector *computeFV(const RootedROMolPtr *ion, const RootedROMolPtr *nl);

  bool includesFeature(const std::string &fname);

private:
  // List of feature classes ready to be used
  static const boost::ptr_vector<Feature> &featureCogs();

  // Indexes of feature classes that are selected for use
  std::vector<int> used_feature_idxs;

  // Helper function - Configure feature for use
  void configureFeature(std::string &name);
};

#endif // __FEATURE_H__