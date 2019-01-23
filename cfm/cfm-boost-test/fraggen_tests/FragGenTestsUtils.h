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

#include "MolData.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/filesystem.hpp>

#ifndef CFM_FRAGGENTESTSUTILS_H
#define CFM_FRAGGENTESTSUTILS_H

FragmentGraph *getTestGraph(std::string smiles_or_inchi,
                            int ionization_mode,
                            bool include_h_loss,
                            bool allow_frag_detours = false);

#endif //CFM_FRAGGENTESTSUTILS_H
