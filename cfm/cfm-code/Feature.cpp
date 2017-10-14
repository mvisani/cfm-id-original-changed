/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features.cpp
#
# Description: 	Code for computing features for fragmentations.
#
#				Assume that we have a config file that lists
the feature
#				vectors to compute (line separated text).
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.
# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include <DataStructs/SparseIntVect.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

#include "Feature.h"
#include "Features/BreakAtomPair.h"
#include "Features/BrokenOrigBondType.h"
#include "Features/ExtraRingFeatures.h"
#include "Features/FeatureHelper.h"
#include "Features/FunctionalGroupFeature.h"
#include "Features/GasteigerCharges.h"
#include "Features/HydrogenMovement.h"
#include "Features/HydrogenRemoval.h"
#include "Features/IonExtraFunctionalGroupFeatures.h"
#include "Features/IonFunctionalGroupFeaturesD2.h"
#include "Features/IonFunctionalGroupFeatures.h"
#include "Features/IonFunctionalGroupRootOnlyFeatures.h"
#include "Features/IonicFeatures.h"
#include "Features/IonNeighbourMMFFAtomType.h"
#include "Features/IonRootAtom.h"
#include "Features/IonRootMMFFAtomType.h"
#include "Features/IonRootPairs.h"
#include "Features/IonRootTriples.h"
#include "Features/IonRootTriplesIncludeBond.h"
#include "Features/NeighbourMMFFFeature.h"
#include "Features/NeighbourOrigBondTypes.h"
#include "Features/NLExtraFunctionalGroupFeatures.h"
#include "Features/NLFunctionalGroupFeaturesD2.h"
#include "Features/NLFunctionalGroupFeatures.h"
#include "Features/NLFunctionalGroupRootOnlyFeatures.h"
#include "Features/NLNeighbourMMFFAtomType.h"
#include "Features/NLRootAtom.h"
#include "Features/NLRootMMFFAtomType.h"
#include "Features/NLRootPairs.h"
#include "Features/NLRootTriples.h"
#include "Features/NLRootTriplesIncludeBond.h"
#include "Features/QuadraticFeatures.h"
#include "Features/RadicalFeatures.h"
#include "Features/RingFeatures.h"
#include "Features/RootAtomFeature.h"
#include "Features/RootPathFeature.h"
#include "Features/IonRootEncodingD3.h"
#include "Features/NLRootEncodingD3.h"
#include "Features/IonRootEncodingD3Short.h"
#include "Features/NLRootEncodingD3Short.h"
#include "Features/IonRootEncodingD4.h"
#include "Features/NLRootEncodingD4.h"
#include "Features/IonRootEncodingD4Long.h"
#include "Features/NLRootEncodingD4Long.h"

const boost::ptr_vector<Feature> &FeatureCalculator::featureCogs() {

  static boost::ptr_vector<Feature> cogs;
  static bool initialised = false;

  if (!initialised) {
    cogs.push_back(new BreakAtomPair());
    cogs.push_back(new BrokenOrigBondType());
    cogs.push_back(new NeighbourOrigBondTypes());
    cogs.push_back(new GasteigerCharges());
    cogs.push_back(new HydrogenMovement());
    cogs.push_back(new HydrogenRemoval());
    cogs.push_back(new IonRootAtom());
    cogs.push_back(new NLRootAtom());
    cogs.push_back(new IonicFeatures());
    cogs.push_back(new IonRootPairs());
    cogs.push_back(new IonRootTriples());
    cogs.push_back(new IonFunctionalGroupFeatures());
    cogs.push_back(new NLFunctionalGroupFeatures());
    cogs.push_back(new IonExtraFunctionalGroupFeatures());
    cogs.push_back(new NLExtraFunctionalGroupFeatures());
    cogs.push_back(new IonFunctionalGroupFeaturesD2());
    cogs.push_back(new NLFunctionalGroupFeaturesD2());
    cogs.push_back(new IonFunctionalGroupRootOnlyFeatures());
    cogs.push_back(new NLFunctionalGroupRootOnlyFeatures());
    cogs.push_back(new NLRootPairs());
    cogs.push_back(new NLRootTriples());
    cogs.push_back(new RadicalFeatures());
    cogs.push_back(new RingFeatures());
    cogs.push_back(new ExtraRingFeatures());
    cogs.push_back(new IonRootMMFFAtomType());
    cogs.push_back(new NLRootMMFFAtomType());
    cogs.push_back(new IonNeighbourMMFFAtomType());
    cogs.push_back(new NLNeighbourMMFFAtomType());
    cogs.push_back(new QuadraticFeatures());
    cogs.push_back(new IonRootTriplesIncludeBond());
    cogs.push_back(new NLRootTriplesIncludeBond());
    cogs.push_back(new NLRootEncodingD3());
    cogs.push_back(new IonRootEncodingD3());
    cogs.push_back(new NLRootEncodingD3Short());
    cogs.push_back(new IonRootEncodingD3Short());
    cogs.push_back(new NLRootEncodingD4());
    cogs.push_back(new IonRootEncodingD4());
    cogs.push_back(new NLRootEncodingD4Long());
    cogs.push_back(new IonRootEncodingD4Long());
    initialised = true;
  }
  return cogs;
}

FeatureCalculator::FeatureCalculator(std::string &config_filename) {

  // Read the config file into a set containing the names of the features
  std::ifstream ifs(config_filename.c_str(), std::ifstream::in);

  if (!ifs.good()) {
    std::cout << "Trouble opening feature config file: " << config_filename
              << std::endl;
  }
  while (ifs.good()) {

    std::string name;
    getline(ifs, name);
    boost::trim(name);
    if (name.length() <= 1)
      break;
    configureFeature(name);
  }

  if (used_feature_idxs.size() == 0) {
    std::cout << "Error reading config file, no features found" << std::endl;
    throw(InvalidConfigException());
  }
}

FeatureCalculator::FeatureCalculator(std::vector<std::string> &feature_list) {

  // Find the relevant feature cog for this name
  std::vector<std::string>::iterator itname = feature_list.begin();
  for (; itname != feature_list.end(); ++itname)
    configureFeature(*itname);

  if (used_feature_idxs.size() == 0)
    std::cout << "Warning: No features found in feature list" << std::endl;
}

std::vector<std::string> FeatureCalculator::getFeatureNames() {

  std::vector<std::string> names;
  std::vector<int>::iterator it = used_feature_idxs.begin();
  for (; it != used_feature_idxs.end(); ++it) {
    const Feature *cog = &(featureCogs()[*it]);
    names.push_back(cog->getName());
  }
  return names;
}

const std::vector<std::string> FeatureCalculator::getValidFeatureNames() {

  static std::vector<std::string> output;
  static bool initialised = false;
  if (initialised)
    return output;

  boost::ptr_vector<Feature>::const_iterator it = featureCogs().begin();
  for (; it != featureCogs().end(); ++it)
    output.push_back(it->getName());

  initialised = true;
  return output;
}

void FeatureCalculator::configureFeature(std::string &name) {

  // Find the relevant feature cog for this name
  boost::ptr_vector<Feature>::const_iterator it = featureCogs().begin();
  bool found = false;
  for (int idx = 0; it != featureCogs().end(); ++it, idx++) {
    if (it->getName() == name) {
      used_feature_idxs.push_back(idx);
      found = true;
      break;
    }
  }
  if (!found) {
    std::cout << "Unrecognised feature: " << name << std::endl;
    throw(InvalidConfigException());
  }
}

unsigned int FeatureCalculator::getNumFeatures() {

  unsigned int count = 1; // Bias
  int quadratic = 0;
  std::vector<int>::iterator it = used_feature_idxs.begin();
  for (; it != used_feature_idxs.end(); ++it) {
    count += featureCogs()[*it].getSize();
    if (featureCogs()[*it].getName() == "QuadraticFeatures")
      quadratic = 1;
  }
  if (quadratic)
    count += (count - 1) * (count - 2) / 2;
  return count;
}

FeatureVector *FeatureCalculator::computeFV(const RootedROMolPtr *ion,
                                            const RootedROMolPtr *nl) {

  FeatureVector *fv = new FeatureVector();

  // Add the Bias Feature
  fv->addFeature(1.0);

  // Compute all other features
  std::vector<int>::iterator it = used_feature_idxs.begin();
  for (; it != used_feature_idxs.end(); ++it) {
    try {
      featureCogs()[*it].compute(*fv, ion, nl);
    } catch (std::exception e) {
      std::cout << "Could not compute " << featureCogs()[*it].getName()
                << std::endl;
      throw FeatureCalculationException("Could not compute " +
                                        featureCogs()[*it].getName());
    }
  }

  return fv;
}

bool FeatureCalculator::includesFeature(const std::string &fname) {
  std::vector<int>::iterator it = used_feature_idxs.begin();
  for (; it != used_feature_idxs.end(); ++it) {
    const Feature *cog = &(featureCogs()[*it]);
    if (cog->getName() == fname)
      return true;
  }
  return false;
}

void FeatureVector::addFeature(double value) {
  if (value != 0.0)
    fv.push_back(fv_idx++);
  else
    fv_idx++;
}

void FeatureVector::addFeatureAtIdx(double value, unsigned int idx) {
  if (fv_idx <= idx)
    fv_idx = idx + 1;
  if (value != 0.0)
    fv.push_back(idx);
}

// Helper functions for multiple features
const std::vector<std::string> &Feature::OKsymbols() {

  static std::vector<std::string> x;
  static bool initialised = false;

  if (!initialised) {
    x.push_back("Br");
    x.push_back("C");
    x.push_back("Cl");
    x.push_back("F");
    x.push_back("I");
    x.push_back("N");
    x.push_back("O");
    x.push_back("P");
    x.push_back("S");
    x.push_back("Se");
    x.push_back("Si");

    initialised = true;
  }
  return x;
}

const std::vector<std::string> &Feature::OKSymbolsLess() {

  static std::vector<std::string> x;
  static bool initialised = false;

  if (!initialised) {
    x.push_back("C");
    x.push_back("N");
    x.push_back("O");
    x.push_back("P");
    x.push_back("S");
    x.push_back("X"); // For all other

    initialised = true;
  }
  return x;
}

void Feature::replaceUncommonWithX(std::string &symbol) const {

  // Replace uncommon symbols with X
  if (symbol != "C" && symbol != "N" && symbol != "O" && symbol != "P" &&
      symbol != "S")
    symbol = "X";
}

int Feature::getSymbolsLessIndex(std::string &symbol) const {
  int index = 0;
  if (symbol == "C")
    index = 0;
  if (symbol == "N")
    index = 1;
  if (symbol == "O")
    index = 2;
  if (symbol == "P")
    index = 3;
  if (symbol == "S")
    index = 4;
  if (symbol == "X")
    index = 5;
}