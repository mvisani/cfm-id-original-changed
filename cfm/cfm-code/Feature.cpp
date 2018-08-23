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

#include "Feature.h"

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

unsigned int Feature::GetSizeOfOKSymbolsLess() const {
    return OKSymbolsLess().size();
}

void Feature::replaceUncommonWithX(std::string &symbol) const {

    // Replace uncommon symbols with X
    const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
    for (auto str : *ok_symbols) {
        if (symbol == str) {
            return;
        }
    }
    symbol = "X";
}

// assume last is always "X"
int Feature::getSymbolsLessIndex(const std::string &symbol) const {
    int index = 0;
    const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
    for (auto str : *ok_symbols) {
        if (symbol == str) {
            break;
        }
        index++;
    }
    return index;
}