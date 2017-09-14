/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# RingFeatures.h
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
#pragma once
#include "../Feature.h"

class RingFeatures : public Feature {
public:
	RingFeatures(){ size = 12; name = "RingFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
private:
	//Helper function - compute the distance between two root 
	//atoms in a molecule (assumes ring break)
	int calcRootDistance(const RootedROMolPtr *mol)  const;
};