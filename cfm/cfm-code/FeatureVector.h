//*#########################################################################
//# Mass Spec Prediction and Identification of Metabolites
//#
//# FeatureVector.h
//#
//# Description: 	Code for computing features for fragmentations.
//#
//# Copyright (c) 2018
//# All rights reserved.
//
//# This file is part of the cfm-id project.
//# The contents are covered by the terms of the GNU Lesser General Public
//# License, which is included in the file license.txt, found at the root
//# of the cfm source tree.
//#########################################################################*/

#pragma once

#include <vector>
#include <iostream>

// Structure to hold a sparse computed feature vector
typedef unsigned int feature_t;

class FeatureVector {
public:
    FeatureVector() { fv_idx = 0; };

    FeatureVector(const FeatureVector &old) {
        fv_idx = old.fv_idx;
        fv = old.fv;
        //std::copy(old.fv.begin(), old.fv.end(), fv.begin());
    };

    void addFeature(double value);

    void addFeatureAtIdx(double value, unsigned int idx);

    void addFeatures(const std::vector<double> &values);

    void addFeatures(const std::vector<int> &values);

    unsigned int getTotalLength() const { return fv_idx; };

    feature_t getFeature(int idx) const { return fv[idx]; };

    std::vector<feature_t>::const_iterator getFeatureBegin() const {
        return fv.begin();
    };

    std::vector<feature_t>::const_iterator getFeatureEnd() const {
        return fv.end();
    };

    unsigned int getNumSetFeatures() const { return fv.size(); };

    void writeDebugInfo(std::ostream &out) const;

    bool equals(const FeatureVector & other_fv) const;

private:
    std::vector<feature_t> fv;
    unsigned int fv_idx;
};
