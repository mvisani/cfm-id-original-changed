/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for features.cpp
#
# Author: Felicity Allen, Fei Wang
# Created: November 2012, 2019
#########################################################################*/
#include <boost/test/unit_test.hpp>

#include "Util.h"
#include "MolData.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <boost/filesystem.hpp>

#ifndef CFM_FEATURESTESTS_H
#define CFM_FEATURESTESTS_H

void initMolProps(romol_ptr_t &mol);

// Return feature RootedROMol
RootedROMolPtr getRootedMolPtr(std::string smiles);

// Retrun feature vector by apply break to input ion and nl str
FeatureVector *getFeatureVector(std::string ion_str, std::string nl_str, std::vector<std::string> & fnames);
// Retrun feature vector by apply break to input ion and nl str
FeatureVector *getFeatureVector(std::string ion_str, std::string nl_str, std::string feature_name);

// Retrun feature vector by apply break to input smiles
FeatureVector *
getFeatureVector(std::string smiles_or_inchi, std::vector<std::string> &fnames, int break_id,
                bool include_h_loss = false,
                int ion_idx = 0);

// Retrun feature vector by apply break to input smiles
FeatureVector *
getFeatureVector(std::string smiles_or_inchi, std::string feature_name, int break_id,
                 bool include_h_loss = false,
                 int ion_idx = 0);

#endif //CFM_FEATURESTESTS_H
