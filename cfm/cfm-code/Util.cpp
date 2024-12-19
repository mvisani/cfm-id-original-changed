/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# util.h
#
# Description: 	Useful functions and definitions
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Util.h"
#include "FunctionalGroups.h"
#include <GraphMol/AtomIterators.h>
#include <GraphMol/FragCatalog/FragCatParams.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/inchi.h>
#include <vector>

double getMassTol(double abs_tol, double ppm_tol, double mass) {
	double mass_tol = (mass / 1000000.0) * ppm_tol;
	if (mass_tol < abs_tol) mass_tol = abs_tol;
	return mass_tol;
}

double getMonoIsotopicMass(const romol_ptr_t &mol) {

	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	double mass              = 0.0;
	unsigned int natoms      = mol->getNumAtoms();
	for (int i = 0; i < natoms; i++) {
		RDKit::Atom *atom  = mol->getAtomWithIdx(i);
		std::string symbol = atom->getSymbol();
		mass += pt->getMostCommonIsotopeMass(symbol);
		mass += atom->getTotalNumHs() * pt->getMostCommonIsotopeMass("H");
	}

	// Adjust the mass by one electron according to the charge
	int charge = RDKit::MolOps::getFormalCharge(*mol.get());
	if (charge == 1) mass -= MASS_ELECTRON;
	if (charge == -1) mass += MASS_ELECTRON;

	return mass;
}

// Helper function to find an atom with the given label
RDKit::Atom *getLabeledAtom(const romol_ptr_t &mol, const char *label) {
	RDKit::ROMol::AtomIterator ai;
	int root = 0;
	for (ai = mol.get()->beginAtoms(); ai != mol.get()->endAtoms(); ++ai) {
		(*ai)->getProp(label, root);
		if (root) break;
	}
	if (root)
		return *ai;
	else
		return nullptr;
}

int moleculeHasSingleRadical(const RDKit::ROMol *romol) {

	unsigned int num_radicals = 0;
	for (RDKit::ROMol::ConstAtomIterator ait = romol->beginAtoms(); ait != romol->endAtoms(); ++ait) {
		int ionic_frag_q;
		(*ait)->getProp("IonicFragmentCharge", ionic_frag_q);
		if (ionic_frag_q != 0) continue; // Don't include radicals on ionic fragments
		num_radicals += (*ait)->getNumRadicalElectrons();
	}
	return (num_radicals == 1);
}

int addIonicChargeLabels(RDKit::ROMol *romol) {

	std::vector<int> mapping;
	unsigned int num_frags = RDKit::MolOps::getMolFrags(*romol, mapping);

	int num_ionic = 0;
	for (auto &ai : romol->atoms()) {
		ai->setProp("IonicFragmentCharge", 0);
		if (num_frags > 1 && ai->getDegree() == 0 && ai->getFormalCharge() != 0) {
			ai->setProp("IonicFragmentCharge", ai->getFormalCharge());
			num_ionic++;
		}
	}
	return num_ionic;
}

void alterNumHs(RDKit::Atom *atom, int H_diff) {
	unsigned int nHs = atom->getTotalNumHs();
	atom->setNoImplicit(true);
	if (atom->getAtomicNum() != 6 && nHs !=4){atom->setNumExplicitHs(nHs + H_diff);}
}

romol_ptr_t createMolPtr(const char *smiles_or_inchi) {
	RDKit::RWMol *rwmol;
	if (std::string(smiles_or_inchi).substr(0, 6) == "InChI=") {
		RDKit::ExtraInchiReturnValues rv;
		rwmol = RDKit::InchiToMol(smiles_or_inchi, rv);
	} else
		rwmol = RDKit::SmilesToMol(smiles_or_inchi);
	auto *mol = static_cast<RDKit::ROMol *>(rwmol);
	addIonicChargeLabels(mol);
	return romol_ptr_t(mol);
}

void softmax(std::vector<double> &weights, std::vector<double> &probs) {
	probs.clear();
	double sum = 0.0;
	for (auto weight : weights) {
		double tmp = std::exp(weight);
		probs.push_back(tmp);
		sum += tmp;
	}
	for (auto &prob : probs) { prob /= sum; }
}

void labelNitroGroup(RDKit::ROMol *mol) {
	// NOTE this is a context specific solution for nitro group single bond oxygen
	auto fparams                      = new RDKit::FragCatParams(PI_BOND_FGRPS_PICKLE);
	const RDKit::MOL_SPTR_VECT &fgrps = fparams->getFuncGroups();
	for (auto &fgrp : fgrps) {
		std::string fg_name;
		fgrp->getProp("_Name", fg_name);

		for (auto ai = mol->beginAtoms(); ai != mol->endAtoms(); ++ai) {
			(*ai)->setProp(fg_name, 0);
			(*ai)->setProp(fg_name + "Charge", 0);
		}
		// The format for each match is (queryAtomIdx, molAtomIdx)
		std::vector<RDKit::MatchVectType> fgp_matches;
		RDKit::SubstructMatch(*mol, *fgrp, fgp_matches);
		for (auto &fgp_match : fgp_matches) {
			for (auto &match : fgp_match) {
				mol->getAtomWithIdx(match.second)->setProp(fg_name, 1);
				mol->getAtomWithIdx(match.second)
				    ->setProp(fg_name + "Charge", mol->getAtomWithIdx(match.second)->getFormalCharge());
			}
		}
	}
	delete fparams;
}

int getValence(const RDKit::Atom *atom) {
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	// Fetch or compute the valence of the atom in the input molecule (we disallow
	// valence changes for now)
	int valence              = -1;
	unsigned int num_val     = pt->getValenceList(atom->getSymbol()).size();
	int def_val              = pt->getDefaultValence(atom->getSymbol());

	// special case for nitrogroup
	int on_nitro_group;
	atom->getProp("NitroGroup", on_nitro_group);
	if (atom->getSymbol() == "O" && atom->getFormalCharge() == -1 && on_nitro_group) {
		valence = 1;
		return valence;
	}
	if (num_val == 1 && def_val != -1) {
		valence = def_val; // Hack to cover many cases - which can otherwise get
		                   // complicated
	} else {
		// This seems to work in most cases....
		valence = atom->getExplicitValence() + atom->getImplicitValence() + atom->getNumRadicalElectrons();
		if (4 - pt->getNouterElecs(atom->getAtomicNum()) > 0) {
			valence += atom->getFormalCharge();
		} else {
			valence -= atom->getFormalCharge();
		}
	}
	return valence;
}
