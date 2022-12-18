/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FeatureHelper.h
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see
param.cpp.
#
# Copyright (c) 2013,2017
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
#pragma once

#include <exception>
#include <string>
#include <iostream>

#include "../FunctionalGroups.h"
#include "../Feature.h"
#include "../FeatureCalculator.h"

class FeatureHelperException : public std::exception {
private:
    std::string message_;

public:
    explicit FeatureHelperException(const std::string &message) noexcept
            : message_(message) {};

    const char *what() const noexcept override {
        std::cout << "Error in FeatureHelper: " << message_ << std::endl;
        return message_.c_str();
    }

    ~FeatureHelperException() noexcept override = default;;
};

class FeatureHelper {
public:
    FeatureHelper() = default;
    explicit FeatureHelper(FeatureCalculator *fc) {
        if(fc != nullptr) {
            exec_flags[0] = fc->includesFeature("GasteigerCharges");
            exec_flags[1] = fc->includesFeature("HydrogenMovement") ||
                            fc->includesFeature("HydrogenRemoval");

            exec_flags[2] = fc->includesFeature("IonRootMMFFAtomType") ||
                            fc->includesFeature("NLRootMMFFAtomType") ||
                            fc->includesFeature("IonNeighbourMMFFAtomType") ||
                            fc->includesFeature("NLNeighbourMMFFAtomType");

            exec_flags[3] = fc->includesFeature("BrokenOrigBondType") ||
                            fc->includesFeature("NeighbourOrigBondTypes") ||
                            fc->includesFeature("BreakHistoryFeature");

            exec_flags[4] = fc->includesFeature("IonFunctionalGroupFeatures") ||
                            fc->includesFeature("NLFunctionalGroupFeatures") ||
                            fc->includesFeature("IonFunctionalGroupFeaturesD2") ||
                            fc->includesFeature("NLFunctionalGroupFeaturesD2") ||
                            fc->includesFeature("IonFunctionalGroupRootOnlyFeatures") ||
                            fc->includesFeature("NLFunctionalGroupRootOnlyFeatures");

            exec_flags[5] = fc->includesFeature("IonExtraFunctionalGroupFeatures") ||
                            fc->includesFeature("NLExtraFunctionalGroupFeatures");

            if (exec_flags[4]) {
                fparams = new RDKit::FragCatParams(FGRPS_PICKLE);
                if (fparams->getNumFuncGroups() != NUM_FGRPS)
                    throw FeatureHelperException(
                            "Mismatch in expected and found number of functional groups");
            }

            if (exec_flags[5]) {
                xfparams = new RDKit::FragCatParams(EXTRA_FGRPS_PICKLE);
                if (xfparams->getNumFuncGroups() != NUM_EXTRA_FGRPS)
                    throw FeatureHelperException(
                            "Mismatch in expected and found number of extra functional groups");
            }
        }
    };

    ~FeatureHelper() {
        delete fparams;
        delete xfparams;
    }

    void addLabels(RDKit::RWMol *rwmol) {
        initialiseRoots(rwmol);
        labelAromatics(rwmol);
        if (exec_flags[0])
            labelGasteigers(rwmol);
        if (exec_flags[1])
            labelOriginalMasses(rwmol);
        if (exec_flags[2])
            labelMMFFAtomTypes(rwmol);
        if (exec_flags[3])
            labelOriginalBondTypes(rwmol);
        if (exec_flags[4])
            labelFunctionalGroups(rwmol, false);
        if (exec_flags[5])
            labelFunctionalGroups(rwmol, true);

        labelAtomsWithLonePairs(rwmol);
    };

    bool getExecFlag(unsigned int idx) { return exec_flags[idx]; };

    static int getBondTypeAsInt(RDKit::Bond *bond);

private:
    std::vector<bool> exec_flags = {false,false,false,false,false,false};

    // Helper functions - used to create labels on atoms and bonds,
    // that will be used in Feature Calculations and can't be computed once
    // a molecule is broken
    static void initialiseRoots(RDKit::RWMol *rwmol);

    static void labelGasteigers(RDKit::RWMol *rwmol);

    static void labelAromatics(RDKit::RWMol *rwmol);

    static void labelOriginalMasses(RDKit::RWMol *rwmol);

    static void labelMMFFAtomTypes(RDKit::RWMol *rwmol);

    static void labelAtomsWithLonePairs(RDKit::RWMol *rwmol);

    static void labelOriginalBondTypes(RDKit::RWMol *rwmol);

    void labelFunctionalGroups(RDKit::RWMol *rwmol,
                               bool extra); // Not static because it uses fparams.

    RDKit::FragCatParams *fparams = nullptr;
    RDKit::FragCatParams *xfparams = nullptr;
};
