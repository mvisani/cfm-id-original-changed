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

#pragma once

#include "Util.h"
#include "FunctionalGroups.h"
#include "FeatureVector.h"

#include <GraphMol/FragCatalog/FragCatParams.h>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/lexical_cast.hpp>

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

typedef std::pair<std::string, std::string> symbol_pair_t;

// Base class to compute a feature - all features should inherit from this

class Feature {

public:
    unsigned int getSize() const { return size; };

    std::string getName() const { return name; };

    virtual ~Feature() {};

protected:
    unsigned int size;
    std::string name;
};

class BreakFeature: public Feature{
public:
    virtual void compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const = 0;

protected:

    static const std::vector<std::string> &OKsymbols();

    static const std::vector<std::string> &OKSymbolsLess();

    void replaceUncommonWithX(std::string &symbol) const;

    int getSymbolsLessIndex(const std::string &symbol) const;

    unsigned int GetSizeOfOKSymbolsLess() const;
};

class FragmentFeature: public Feature{
public:
    virtual void compute(FeatureVector &fv, romol_ptr_t precursor_ion) const = 0;
};