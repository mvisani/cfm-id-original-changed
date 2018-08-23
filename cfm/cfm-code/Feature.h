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
    virtual void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const = 0;

    unsigned int getSize() const { return size; };

    std::string getName() const { return name; };

    virtual ~Feature() {};

protected:
    unsigned int size;
    std::string name;

    static const std::vector<std::string> &OKsymbols();

    static const std::vector<std::string> &OKSymbolsLess();

    void replaceUncommonWithX(std::string &symbol) const;

    int getSymbolsLessIndex(const std::string &symbol) const;

    unsigned int GetSizeOfOKSymbolsLess() const;
};

/*class BreakFeature: public Feature{
    virtual void compute(FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl, int depth) const = 0;
};

class FragmentFeature: public Feature{
    virtual void compute(FeatureVector &fv, std::string fragment_smiles) const = 0;
};*/