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
#include "FeaturesTestsUtil.h"

void initMolProps(romol_ptr_t &mol) {
    RDKit::ROMol::AtomIterator ai;
    for (ai = mol.get()->beginAtoms(); ai != mol.get()->endAtoms(); ++ai) {
        (*ai)->setProp("Root", 0);
    }
    mol.get()->setProp("IsRingBreak", 0);
}

RootedROMol getRootedMolPtr(std::string smiles) {
    romol_ptr_t mol = createMolPtr(smiles.c_str());
    initMolProps(mol);
    RootedROMol rtd_mol(mol, mol.get()->getAtomWithIdx(0));
    return rtd_mol;
}

// Retrun feature vector by apply break to input ion and nl str
FeatureVector *getFeatureVector(std::string ion_str, std::string nl_str, std::vector<std::string> & fnames) {

    FeatureCalculator *fc = new FeatureCalculator(fnames);

    // init ion
    romol_ptr_t ion = createMolPtr(ion_str.c_str());
    initMolProps(ion);
    RootedROMol rtd_ion(ion, ion.get()->getAtomWithIdx(0));

    // init nl
    romol_ptr_t nl = createMolPtr(nl_str.c_str());
    initMolProps(nl);
    RootedROMol rtd_nl(nl, nl.get()->getAtomWithIdx(0));

    FeatureVector *fv = fc->computeFeatureVector(&rtd_ion, &rtd_nl, nullptr);

    delete fc;
    return fv;
}

// Retrun feature vector by apply break to input ion and nl str
FeatureVector *getFeatureVector(std::string ion_str, std::string nl_str, std::string feature_name) {

    std::vector<std::string> fnames;
    fnames.push_back(feature_name);

    return getFeatureVector(ion_str, nl_str, fnames);
}

// Retrun feature vector by apply break to input smiles
FeatureVector *
getFeatureVector(std::string smiles_or_inchi, std::vector<std::string> &fnames, int break_id,
                bool include_h_loss, int ion_idx) {

    FeatureCalculator *fc = new FeatureCalculator(fnames);

    FragmentGraphGenerator fgen(fc);
    FragmentTreeNode *node = fgen.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
    std::vector<Break> breaks;
    node->generateBreaks(breaks, include_h_loss, false);

    //std::cout << breaks.size() << std::endl;

    node->applyBreak(breaks[break_id], ion_idx);
    node->generateChildrenOfBreak(breaks[break_id], false);
    FragmentTreeNode *child = &(node->children[0]);
    Transition transition(-1, -1, child->nl, child->ion);

    FeatureVector *fv = fc->computeFeatureVector(transition.getIon(), transition.getNeutralLoss(), nullptr);

    delete fc;
    return fv;
}

// Retrun feature vector by apply break to input smiles
FeatureVector *
getFeatureVector(std::string smiles_or_inchi, std::string feature_name, int break_id, bool include_h_loss,
        int ion_idx) {

    std::vector<std::string> fnames;
    fnames.push_back(feature_name);

    return getFeatureVector(smiles_or_inchi, fnames, break_id, include_h_loss, ion_idx);
}