/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# BreakAtomPair.cpp
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see param.cpp.
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#include "BreakAtomPair.h"

void


BreakAtomPair::compute(FeatureVector &fv, const RootedROMol *ion, const RootedROMol *nl) const {

    int ring_break;
    nl->mol.get()->getProp("IsRingBreak", ring_break);

    //Ion Symbol
    std::string ion_root_symbol = ion->root->getSymbol();
    replaceUncommonWithX(ion_root_symbol, false);

    //Neutral Loss Symbol
    std::string nl_root_symbol = nl->root->getSymbol();
    replaceUncommonWithX(nl_root_symbol, false);

    //Pairs
    symbol_pair_t pair(ion_root_symbol, nl_root_symbol);

    //Iterate through all combinations of atom pairs, appending
    //a feature for each; 1 if it matches, 0 otherwise.
    //Note: the order matters here, ion first then nl
    std::vector<std::string>::const_iterator it1, it2;
    const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
    for (it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1) {
        for (it2 = ok_symbols->begin(); it2 != ok_symbols->end(); ++it2) {
            symbol_pair_t sp = symbol_pair_t(*it1, *it2);
            double feature = 0.0;
            if (sp == pair)
                feature = 1.0;
            fv.addFeature(feature);
        }
    }

    if(ring_break)
        fv.addFeature(1);
    else
        fv.addFeature(0);
}
