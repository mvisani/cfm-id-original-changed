/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraph.cpp
#
# Description: 	FragmentGraph class for holding the results of a generated
#				fragment graph.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FragmentGraph.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/MolOps.h>
#include <stack>

Transition::Transition(int a_from_id, int a_to_id, const romol_ptr_t &a_nl, const romol_ptr_t &an_ion) {

    RDKit::Atom *root = nullptr;
    RDKit::Atom *first_atom = an_ion.get()->getAtomWithIdx(0);
    if (first_atom->hasProp("Root"))
        root = getLabeledAtom(an_ion, "Root");
    else
        std::cerr << "Warning: Ion does not have Root Prop" << std::endl;
    if (root == nullptr)
       std::cerr << "Warning: Ion Root atoms not defined" << std::endl;
 
    ion = RootedROMol(an_ion, root);

    root = nullptr;
    first_atom = a_nl.get()->getAtomWithIdx(0);
    if (first_atom->hasProp("Root"))
        root = getLabeledAtom(a_nl, "Root");
    else
       std::cerr << "Warning: NL does not have Root Prop" << std::endl;
    if (root == nullptr)
       std::cerr << "Warning: NL Root atoms not defined" << std::endl;
    nl = RootedROMol(a_nl, root);

    to_id = a_to_id;
    from_id = a_from_id;
    nl_smiles = RDKit::MolToSmiles(*(nl.mol.get()));

}

Transition::Transition(int a_from_id, int a_to_id, const RootedROMol &a_nl, const RootedROMol &an_ion)
        :
        to_id(a_to_id), from_id(a_from_id), nl(a_nl), ion(an_ion) {
    nl_smiles = RDKit::MolToSmiles(*(nl.mol.get()));

}

//Add a fragment node to the graph (should be the only way to modify the graph)
//	-- Add a fragment, or return an id, if it already exists
//  -- Add a transition to this fragment, based on the provided parent fragment id
//  -- Update the relevant tmaps
// Note: If parentid < 0 (i.e. starting node) doesn't add transition.
int FragmentGraph::addToGraph(const FragmentTreeNode &node, int parentid) {

    //If the fragment doesn't exist, add it
    double mass = getMonoIsotopicMass(node.ion);
    int id = addFragmentOrFetchExistingId(node.ion, mass, node.isIntermediate(), node.isCyclization());

    if (parentid < 0 || fragments[id]->getDepth() == -1)
        fragments[id]->setDepth(node.depth);    //Set start fragment depth

    //Add a transition, IF:
    // - This transition does not exist (i.e. this parent_id -> id)
    // - There is no (already added) shorter path to this fragment (we track the minimum depth at which each fragment is seen)
    //  (NOTE: This will not remove previously added longer paths, we will remove them all at once at the end)
    if (parentid >= 0 && findMatchingTransition(parentid, id) < 0 &&
        (allow_frag_detours || node.depth <= fragments[id]->getDepth())) {

        if (node.depth < fragments[id]->getDepth())
            fragments[id]->setDepth(node.depth);    //Update the depth of the fragment

        int idx = (int) transitions.size();
        transitions.push_back(std::make_shared<Transition>(parentid, id, node.nl, node.ion));

        //Update the tmaps
        from_id_tmap[parentid].push_back(idx);
        to_id_tmap[id].push_back(idx);
    }

    return id;
}

//As for previous function, but delete the mols in the transition and compute and store a feature vector instead
int FragmentGraph::addToGraphAndReplaceMolWithFV(const FragmentTreeNode &node, int parent_frag_id, FeatureCalculator *fc) {

    //If the fragment doesn't exist, add it
    double mass = getMonoIsotopicMass(node.ion);
    int frag_id = addFragmentOrFetchExistingId(node.ion, mass, node.isIntermediate(), node.isCyclization());

    if (fragments[frag_id]->getDepth() == -1 || node.depth < fragments[frag_id]->getDepth())
        fragments[frag_id]->setDepth(node.depth);    //Set start fragment depth

    // find if this trans does exist
    int existing_trans_id = -1;
    if (parent_frag_id >= 0)
        existing_trans_id = findMatchingTransition(parent_frag_id, frag_id);

    //Add a transition, if one does not exist
    if (parent_frag_id >= 0 && existing_trans_id < 0 &&
        (allow_frag_detours || node.depth <= fragments[frag_id]->getDepth())) {

        int idx = transitions.size();
        transitions.push_back(std::make_shared<Transition>(parent_frag_id, frag_id, node.nl, node.ion));

        //Update the tmaps
        from_id_tmap[parent_frag_id].push_back(idx);
        to_id_tmap[frag_id].push_back(idx);

        //Compute a feature vector
        auto t = transitions.back();
        FeatureVector *fv;
        std::string frag_smiles = *fragments[parent_frag_id]->getIonSmiles();
        auto frag_ptr = createMolPtr(frag_smiles.c_str());

        try {
            fv = fc->computeFeatureVector(t->getIon(), t->getNeutralLoss(), frag_ptr);
            t->setFeatureVector(fv);
        }
        catch (FeatureCalculationException &e) {
            //If we couldn't compute the feature vector, set a dummy feature vector with bias only.
            fv = new FeatureVector();
            fv->addFeatureAtIdx(1.0, idx);
            fv->addFeatureAtIdx(0.0, fc->getNumFeatures() - 1);
        }

        //Delete the mols
        t->deleteIon();
        t->deleteNeutralLoss();
    } else if (parent_frag_id >= 0 && existing_trans_id >= 0 &&
               (node.depth == fragments[frag_id]->getDepth())) {
        // if those transitions does exists, create a duplication
        int trans_idx = transitions.size();
        transitions.push_back(std::make_shared<Transition>());
        auto trans = transitions.back();
        trans->createdDuplication(*transitions[existing_trans_id]);
        //Update the tmaps
        from_id_tmap[parent_frag_id].push_back(trans_idx);
        to_id_tmap[frag_id].push_back(trans_idx);
    }

    return frag_id;
}

int
FragmentGraph::addToGraphWithThetas(const FragmentTreeNode &node, const std::vector<double> *thetas,
                                    int parent_frag_id) {

    //If the fragment doesn't exist, add it
    double mass = getMonoIsotopicMass(node.ion);
    int frag_id = addFragmentOrFetchExistingId(node.ion, mass, node.isIntermediate(), node.isCyclization());

    if (parent_frag_id < 0 || fragments[frag_id]->getDepth() == -1)
        fragments[frag_id]->setDepth(node.depth);    //Set start fragment depth

    //Add a transition, if one does not exist
    int existing_trans_id = -1;
    if (parent_frag_id >= 0)
        existing_trans_id = findMatchingTransition(parent_frag_id, frag_id);

    if (parent_frag_id >= 0 && existing_trans_id < 0 &&
        (allow_frag_detours || node.depth <= fragments[frag_id]->getDepth())) {

        int idx = transitions.size();
        transitions.push_back(std::make_shared<Transition>(parent_frag_id, frag_id, node.nl, node.ion));

        //Update the tmaps
        from_id_tmap[parent_frag_id].push_back(idx);
        to_id_tmap[frag_id].push_back(idx);

        //Set the theta values and delete the mols
        auto t = transitions.back();
        t->setTmpThetas(thetas);
        t->deleteIon();
        t->deleteNeutralLoss();

    } else if (parent_frag_id >= 0 && existing_trans_id >= 0 &&
               (node.depth == fragments[frag_id]->getDepth())) {
        // if those transitions does exist, create a duplication
        int trans_idx = transitions.size();
        transitions.push_back(std::make_shared<Transition>());
        auto trans = transitions.back();
        trans->createdDuplication(*transitions[existing_trans_id]);
        //Update the tmaps
        from_id_tmap[parent_frag_id].push_back(trans_idx);
        to_id_tmap[frag_id].push_back(trans_idx);
    }

    return frag_id;
}


int FragmentGraph::addFragmentOrFetchExistingId(romol_ptr_t ion, double mass, 
    bool is_intermediate,
    bool is_cyclization) {

    std::string reduced_smiles;
    //Round the mass to 5 decimal places and use that as an initial filter
    double rounded_mass = floor(mass * 10000.0 + 0.5) / 10000.0;
    if (frag_mass_lookup.find(rounded_mass) != frag_mass_lookup.end()) {

        //Create a copy of the ion and then reduce it, making all bonds single and filling in hydrogens
        RDKit::RWMol f1_copy = *ion.get();
        reduceMol(f1_copy);
        RDKit::MolOps::sanitizeMol(f1_copy);

        //Found an entry with this mass, check the linked fragments for a match
        std::vector<int>::iterator it;
        it = frag_mass_lookup[rounded_mass].begin();
        for (; it != frag_mass_lookup[rounded_mass].end(); ++it) {
            try {
                RDKit::RWMol *f2_reduced = RDKit::SmilesToMol(*fragments[*it]->getReducedSmiles());

                if (areMatching(&f1_copy, f2_reduced)) {
                    delete f2_reduced;
                    return *it;
                }
                delete f2_reduced;
            }
            catch (RDKit::MolSanitizeException &e) {
                std::cout << "Could not sanitize " << *fragments[*it]->getReducedSmiles() << std::endl;
                throw &e;
            }
        }
        reduced_smiles = RDKit::MolToSmiles(f1_copy);
    } else
        frag_mass_lookup[rounded_mass];

    //No match found, create the fragment
    if (reduced_smiles.size() == 0) {
        RDKit::RWMol f1_copy = *ion.get();
        reduceMol(f1_copy);
        reduced_smiles = RDKit::MolToSmiles(f1_copy);
    }
    std::string smiles = RDKit::MolToSmiles(*ion.get());
    int newid = fragments.size();
    if (include_isotopes) {
        Spectrum isotope_spectrum;
        long charge = RDKit::MolOps::getFormalCharge(*ion.get());
        isotope->computeIsotopeSpectrum(isotope_spectrum, ion, charge);
        fragments.push_back(new Fragment(smiles, reduced_smiles, newid, mass, isotope_spectrum, is_intermediate, is_cyclization));
    } else{
        fragments.push_back(new Fragment(smiles, reduced_smiles, newid, mass, is_intermediate, is_cyclization));
    }
        

    frag_mass_lookup[rounded_mass].push_back(newid);
    from_id_tmap.resize(newid + 1);
    to_id_tmap.resize(newid + 1);
    return newid;
}

bool FragmentGraph::areMatching(RDKit::ROMol *f1_reduced_ion, RDKit::ROMol *f2_reduced_ion) {

    //Quick preliminary check to throw away non-matches
    if (f1_reduced_ion->getNumAtoms() != f2_reduced_ion->getNumAtoms())
        return false;

    //Note: this will fail if one mol is a substruct
    //of another but not an exact match, however given we've just
    //checked that the number of atoms match, this can't happen
    RDKit::MatchVectType match;
    return RDKit::SubstructMatch(*f1_reduced_ion, *f2_reduced_ion, match, false, false);
}

void FragmentGraph::reduceMol(RDKit::RWMol &rwmol) {

    //Ensure that the OrigValence tags are set, otherwise add them
    //int israd = moleculeHasSingleRadical(&rwmol);
    /*RDKit::ROMol::AtomIterator ai;
    if (!rwmol.getAtomWithIdx(0)->hasProp("OrigValence")) {
        std::cout << "OrigValence not set" << std::endl;
        for (ai = rwmol.beginAtoms(); ai != rwmol.endAtoms(); ++ai) {
            int origval = (*ai)->getExplicitValence() + (*ai)->getImplicitValence() + (*ai)->getNumRadicalElectrons() -
                          (*ai)->getFormalCharge();

            // For radicals generated by removing a multibond,
            // charge and radical are separated,
            // so charge needs to be added, not subtracted
            // sometime this fails ... we need to check if it is multibond
            if (israd
                && (*ai)->getFormalCharge()
                && !(*ai)->getNumRadicalElectrons()
                && (*ai)->getSymbol() == "C"){
                origval += 2 * (*ai)->getFormalCharge();
            }

            (*ai)->setProp("OrigValence", origval);
        }
    }*/

    //Set all bonds to SINGLE
    RDKit::ROMol::BondIterator bi;
    for (bi = rwmol.beginBonds(); bi != rwmol.endBonds(); ++bi) {
        (*bi)->setBondType(RDKit::Bond::SINGLE);
        (*bi)->setIsAromatic(false);
    }

    //Set all Hydrogens to give full valence and no charge or radicals
    RDKit::ROMol::AtomIterator ai;
    for (ai = rwmol.beginAtoms(); ai != rwmol.endAtoms(); ++ai) {
        (*ai)->setFormalCharge(0);
        (*ai)->setNumRadicalElectrons(0);
        int valence;
        (*ai)->getProp("OrigValence", valence);
        int numH = valence - (*ai)->getDegree();
        if (numH < 0) {
            (*ai)->setFormalCharge(1);
            numH = 0;
        }
        (*ai)->setNumExplicitHs(numH);
        (*ai)->setIsAromatic(false);
    }

}

//Find the id for an existing transition that matches the input ids
//or -1 in the case where no such transition is found
int FragmentGraph::findMatchingTransition(int from_id, int to_id) {
    std::vector<int>::iterator it = from_id_tmap[from_id].begin();
    for (; it != from_id_tmap[from_id].end(); ++it) {
        if (*it >= transitions.size())
            continue;
        if (transitions[*it]->getToId() == to_id)
            return *it;
    }
    return -1;
}

//Write the Fragments only to file (formerly the backtrack output - without extra details)
void FragmentGraph::writeFragmentsOnly(std::ostream &out) const {

    for (auto &it : fragments) {
        out << it->getId() << " ";
        out << std::setprecision(6) << it->getMass() << " ";
        out << *(it->getIonSmiles());

        if (it->isIntermediate())
            out << " Intermediate Fragment";
        if (it->isCyclization())
            out << " Cyclization Fragment";
        out << std::endl;
    }
}

void FragmentGraph::writeFragmentsOnlyForIds(std::ostream &out, std::set<int> & ids) const {

    for (auto &it : fragments) {
        if (ids.find(it->getId()) != ids.end()) {
            out << it->getId() << " ";
            out << std::setprecision(6) << it->getMass() << " ";
            out << *(it->getIonSmiles()) << std::endl;
        }
    }


}

//Write the FragmentGraph to file (formerly the transition output - without feature details)
void FragmentGraph::writeFullGraph(std::ostream &out) const {

    //Fragments
    out << fragments.size() << std::endl;
    writeFragmentsOnly(out);
    out << std::endl;

    //Transitions
    for (const auto &transition : transitions) {
        if (!transition->isDuplicate()) {
            out << transition->getFromId() << " ";
            out << transition->getToId() << " ";
            out << *fragments[transition->getToId()]->getIonSmiles() << " ";
            out << *(transition->getNLSmiles()) << " ";
            out << std::endl;
        }
    }
}

void FragmentGraph::writeFeatureVectorGraph(std::ostream &out, bool include_isotopes) const {

    // Write Depth of Graph
    out.write(reinterpret_cast<const char *>(&depth), sizeof(depth));
    //Fragments
    unsigned int inlcis = include_isotopes;
    out.write(reinterpret_cast<const char *>(&inlcis), sizeof(inlcis));
    unsigned int numf = fragments.size();
    out.write(reinterpret_cast<const char *>(&numf), sizeof(numf));
    auto it = fragments.begin();
    for (; it != fragments.end(); ++it) {
        // write id
        int id = (*it)->getId();
        out.write(reinterpret_cast<const char *>(&id), sizeof(id));
        // write mass
        double mass = (*it)->getMass();
        out.write(reinterpret_cast<const char *>(&mass), sizeof(mass));

        // write flag
        bool is_intermediate = (*it)->isIntermediate();
        out.write(reinterpret_cast<const char *>(&is_intermediate), sizeof(is_intermediate));

        // write flag
        bool is_cyclization = (*it)->isCyclization();
        out.write(reinterpret_cast<const char *>(&is_cyclization), sizeof(is_cyclization));

        if (include_isotopes) {
            const Spectrum *isospec = (*it)->getIsotopeSpectrum();
            unsigned int iso_size = isospec->size();
            out.write(reinterpret_cast<const char *>(&iso_size), sizeof(iso_size));
            Spectrum::const_iterator itp = isospec->begin();
            for (; itp != isospec->end(); ++itp) {
                out.write(reinterpret_cast<const char *>(&itp->mass), sizeof(itp->mass));
                out.write(reinterpret_cast<const char *>(&itp->intensity), sizeof(itp->intensity));
            }
        }
    }

    //Transitions
    unsigned int numt = transitions.size();
    out.write(reinterpret_cast<const char *>(&numt), sizeof(numt));
    auto itt = transitions.begin();
    for (; itt != transitions.end(); ++itt) {
        int fromid = (*itt)->getFromId();
        int toid = (*itt)->getToId();
        out.write(reinterpret_cast<const char *>(&fromid), sizeof(fromid));
        out.write(reinterpret_cast<const char *>(&toid), sizeof(toid));
        FeatureVector *fv = (*itt)->getFeatureVector();
        unsigned int num_set = fv->getNumSetFeatures();
        unsigned int fv_len = fv->getTotalLength();
        out.write(reinterpret_cast<const char *>(&num_set), sizeof(num_set));
        out.write(reinterpret_cast<const char *>(&fv_len), sizeof(fv_len));

        for (auto fv_it = fv->getFeatureBegin(); fv_it != fv->getFeatureEnd(); ++fv_it) {
            unsigned int fv_idx = *fv_it;
            out.write(reinterpret_cast<const char *>(&fv_idx), sizeof(fv_idx));
        }
    }
}

void FragmentGraph::readFeatureVectorGraph(std::istream &ifs) {

    std::string null = "";
    // Set Depth of Graph
    ifs.read(reinterpret_cast<char *>(&depth), sizeof(depth));
    // Include Isotope Flag
    unsigned int include_isotopes;
    ifs.read(reinterpret_cast<char *>(&include_isotopes), sizeof(include_isotopes));
    // Fragaments
    unsigned int numf;
    ifs.read(reinterpret_cast<char *>(&numf), sizeof(numf));
    for (int i = 0; i < numf; i++) {

        int id;
        ifs.read(reinterpret_cast<char *>(&id), sizeof(id));
        double mass;
        ifs.read(reinterpret_cast<char *>(&mass), sizeof(mass));

        bool is_intermediate;
        ifs.read(reinterpret_cast<char *>(&is_intermediate), sizeof(is_intermediate));

        bool is_cyclization;
        ifs.read(reinterpret_cast<char *>(&is_cyclization), sizeof(is_cyclization));

        if (include_isotopes) {
            Spectrum isospec;
            unsigned int iso_size;
            ifs.read(reinterpret_cast<char *>(&iso_size), sizeof(iso_size));
            for (int j = 0; j < iso_size; j++) {
                double pmass;
                ifs.read(reinterpret_cast<char *>(&pmass), sizeof(pmass));
                double pintensity;
                ifs.read(reinterpret_cast<char *>(&pintensity), sizeof(pintensity));
                isospec.push_back(Peak(pmass, pintensity));
            }
            fragments.push_back(new Fragment(null, null, id, mass, isospec, is_intermediate, is_cyclization));
        } else
            fragments.push_back(new Fragment(null, null, id, mass, is_intermediate, is_cyclization));
    }

    //Transitions
    unsigned int numt;
    ifs.read(reinterpret_cast<char *>(&numt), sizeof(numt));
    for (int i = 0; i < numt; i++) {

        int fromid;
        ifs.read(reinterpret_cast<char *>(&fromid), sizeof(fromid));
        int toid;
        ifs.read(reinterpret_cast<char *>(&toid), sizeof(toid));

        FeatureVector *fv = new FeatureVector();
        unsigned int num_set;
        ifs.read(reinterpret_cast<char *>(&num_set), sizeof(num_set));
        unsigned int fv_len;
        ifs.read(reinterpret_cast<char *>(&fv_len), sizeof(fv_len));
        for (unsigned int j = 0; j < num_set; j++) {
            unsigned int fidx;
            ifs.read(reinterpret_cast<char *>(&fidx), sizeof(fidx));
            fv->addFeatureAtIdx(1.0, fidx);
        }
        if (fv->getTotalLength() != fv_len)
            fv->addFeatureAtIdx(0.0, fv_len - 1);

        transitions.push_back(std::make_shared<Transition>(fromid, toid, &null));
        transitions[i]->setFeatureVector(fv);
    }

    //Create from_id and to_id maps
    to_id_tmap.resize(fragments.size());
    from_id_tmap.resize(fragments.size());
    auto it = transitions.begin();
    for (int idx = 0; it != transitions.end(); ++it, idx++) {
        to_id_tmap[(*it)->getToId()].push_back(idx);
        from_id_tmap[(*it)->getFromId()].push_back(idx);
    }
}

//Function to remove detour transitions from the graph (used if !cfg.allow_frag_detours)
void FragmentGraph::removeDetours() {
    std::vector<int> remove_ids;
    // Get a list of  transitions need to be removed
    for (int i = 0; i < transitions.size(); i++) {
        int from_id = transitions[i]->getFromId();
        int to_id = transitions[i]->getToId();
        if (fragments[from_id]->getDepth() >= fragments[to_id]->getDepth())
            remove_ids.push_back(i);
    }

    // Copy and pasted code
    // this is urgly but seems there is better way
    std::sort(remove_ids.begin(), remove_ids.end());
    //std::cout << transitions.size() << std::endl;
    std::vector<int> id_map(transitions.size());    //Map old id to new id (or -1 if deleting)
    //Remove transitions, and record a mapping of old->new transition ids
    unsigned int next_id = 0;
    unsigned int id_idx = 0;
    for (int i = 0; i < transitions.size(); i++) {
        bool need_remove = (id_idx < remove_ids.size());
        if (need_remove)
            need_remove = (need_remove && (remove_ids[id_idx] == i));
        // if deleting
        if (need_remove) {
            id_map[i] = -1;
            ++id_idx;
        } else {
            transitions[next_id] = transitions[i];
            id_map[i] = next_id++;
        }
    }
    transitions.resize(next_id);

    //Update the id maps
    for (int i = 0; i < fragments.size(); i++) {

        //Update the To Id Tmap
        unsigned int count = 0;
        for (int j = 0; j < to_id_tmap[i].size(); j++) {
            if (id_map[to_id_tmap[i][j]] >= 0)
                to_id_tmap[i][count++] = id_map[to_id_tmap[i][j]];
        }
        to_id_tmap[i].resize(count);

        //Update the From Id Tmap
        count = 0;
        for (int j = 0; j < from_id_tmap[i].size(); j++) {
            if (id_map[from_id_tmap[i][j]] >= 0)
                from_id_tmap[i][count++] = id_map[from_id_tmap[i][j]];
        }
        from_id_tmap[i].resize(count);
    }
}

void
FragmentGraph::getSampledTransitionIdsRandomWalk(std::set<int> &selected_ids, int max_selection) {

    // use ceil so we are at least get one

    std::vector<std::uniform_int_distribution<int>> uniform_int_distributions;
    // add to discrete_distributions
    for (auto &frag_trans_ids: from_id_tmap)
        uniform_int_distributions.emplace_back(std::uniform_int_distribution<>(0, (int) frag_trans_ids.size() - 1));

    for(int i = 0; i < max_selection; ++i) {
        // init queue and add root
        std::queue<int> fgs;
        std::set<int> visited_fgs;
        fgs.push(0);

        while (!fgs.empty()) {
            // get current id
            int frag_id = fgs.front();
            visited_fgs.insert(frag_id);
            fgs.pop();

            // if there is somewhere to go
            if (!from_id_tmap[frag_id].empty()) {
                // add a uct style random select
                int selected_idx = uniform_int_distributions[frag_id](util_rng);
                if (selected_idx < from_id_tmap[frag_id].size()) {
                    // go to child
                    int selected_trans_id = from_id_tmap[frag_id][selected_idx];
                    int next_fg_id = transitions[selected_trans_id]->getToId();
                    //make sure we are not visit the same place twice
                    //without this check this may ends in endless loop
                    if (visited_fgs.count(next_fg_id) == 0) {
                        fgs.push(next_fg_id);
                        selected_ids.insert(selected_trans_id);
                    }
                }
            }
        }
    }
}

void FragmentGraph::getSampledTransitionIdsWeightedRandomWalk(std::set<int> &selected_ids,
                                                              int max_num_iter,
                                                              std::vector<double> &thetas,
                                                              double explore_weight) {

    std::vector<std::uniform_int_distribution<int>> uniform_int_distributions;
    for (auto &frag_trans_ids: from_id_tmap)
        uniform_int_distributions.emplace_back(std::uniform_int_distribution<>(0, (int) frag_trans_ids.size() - 1));


    std::vector<std::discrete_distribution<int>> discrete_distributions;

    for (auto &frag_trans_ids: from_id_tmap) {

        // Init weights and prob vector
        std::vector<double> probs;
        std::vector<double> weights;

        for (auto &trans_id : frag_trans_ids) {
            weights.push_back(thetas[trans_id]);
        }

        // Append 0.0 for i -> i
        // exp(0.0) = 1.0
        weights.push_back(0.0);

        // Apply softmax
        softmax(weights, probs);

        // add to discrete_distributions
        discrete_distributions.emplace_back(std::discrete_distribution<int>(probs.begin(), probs.end()));
    }

    int num_iter = 0;
    std::discrete_distribution<int> explore_coin({1.0 - explore_weight, explore_weight});

    while (num_iter <= max_num_iter) {
        // init queue and add root
        std::queue<int> fgs;
        fgs.push(0);

        while (!fgs.empty()) {

            // get current id
            int frag_id = fgs.front();
            fgs.pop();

            // if there is somewhere to go
            if (!from_id_tmap[frag_id].empty()) {
                // add a uct style random select
                int selected_idx = -1;
                int coin = explore_coin(util_rng);
                if (coin == 1) {
                    selected_idx = uniform_int_distributions[frag_id](util_rng);
                } else {
                    selected_idx = discrete_distributions[frag_id](util_rng);
                }
                if (selected_idx < from_id_tmap[frag_id].size() && selected_idx > -1) {
                    // go to child
                    int selected_trans_id = from_id_tmap[frag_id][selected_idx];
                    fgs.push(transitions[selected_trans_id]->getToId());
                    selected_ids.insert(selected_trans_id);
                }
            }
        }
        num_iter++;
    }
}

void FragmentGraph::
getSampledTransitionIdsDiffMapChildOnly(std::set<int> &selected_ids, std::set<unsigned int> &selected_weights){
    std::set<int> visited;
    std::vector<int> path;
    getSampledTransitionIdsWeightDiffChildOnly(selected_weights, visited, 0, path,
                                               selected_ids);
}

void FragmentGraph::
getSampledTransitionIdsDiffMap(std::set<int> &selected_ids, std::set<unsigned int> &selected_weights){
    std::set<int> visited;
    std::map<double, std::set<int>> selected_trans_map;

    std::map<int, std::vector<int>> frag_trans_map;
    std::vector<std::pair<int,int>> frag_trans_pair_path;
    getSampledTransitionIdsWeightDiffs(selected_weights, visited, 0, frag_trans_pair_path, frag_trans_map);

    for(const auto & record : frag_trans_map){
        if(record.second.size() > 1)
            for(const auto & trans_id :record.second)
                selected_ids.insert(trans_id);
    }
}

void FragmentGraph::
getSampledTransitionIdsWeightDiffs(std::set<unsigned int> &selected_weights, std::set<int> &visited,
                                   int frag_id, std::vector<std::pair<int, int>> &path,
                                   std::map<int, std::vector<int>> &trans_to_interest_frags_map) {

    double frag_mass = fragments[frag_id]->getMass();

    if (is_match(selected_weights, frag_mass)) {
        // record troubled frag_id to each transition
        for (const auto &frag_trans_pair : path) {
            if (trans_to_interest_frags_map.find(frag_trans_pair.first) == trans_to_interest_frags_map.end())
                trans_to_interest_frags_map[frag_trans_pair.first];
            trans_to_interest_frags_map[frag_trans_pair.first].push_back(frag_trans_pair.second);
        }
    }

    // if we have see this before
    if (visited.find(frag_id) != visited.end())
        return;

    visited.insert(frag_id);
    for (const auto &trans_id : from_id_tmap[frag_id]) {
        auto current_path = path;
        current_path.push_back(std::pair<int, int>(frag_id, trans_id));
        getSampledTransitionIdsWeightDiffs(selected_weights, visited, transitions[trans_id]->getToId(),
                                           current_path, trans_to_interest_frags_map);
    }

}

void FragmentGraph::
getSampledTransitionIdsWeightDiffChildOnly(std::set<unsigned int> &selected_weights, std::set<int> &visited,
                                           int frag_id, std::vector<int> &path, std::set<int> &selected_ids) {
    // Note since we always start from root
    // and there is no arc lead to root ( or we does not care )
    // we should check all the child fragments

    // if we have see this before
    if (visited.find(frag_id) != visited.end())
        return;

    visited.insert(frag_id);
    // matched ids in childs
    std::vector<int> matched_selected_ids;
    std::vector<int> trans_ids;

    for (const auto &trans_id : from_id_tmap[frag_id]) {
        std::vector<int> path_to_current_child = path;
        path_to_current_child.push_back(trans_id);
        auto child_frag_id = transitions[trans_id]->getToId();
        // check child fragmentation weights
        double child_frag_mass = fragments[child_frag_id]->getMass();

        // check if child mass matches what we are looking for
        if (is_match(selected_weights, child_frag_mass))
            matched_selected_ids.push_back(trans_id);

        getSampledTransitionIdsWeightDiffChildOnly(selected_weights, visited, child_frag_id,
                                                   path_to_current_child,
                                                   selected_ids);
    }

    if (!matched_selected_ids.empty()) {

        // add all trans lead to this node
        for (const auto &trans_id : path)
            selected_ids.insert(trans_id);

        //add all trans in this node
        for (const auto &trans_id : matched_selected_ids)
           selected_ids.insert(trans_id);
    }
}

bool FragmentGraph::is_match(std::set<unsigned int> &weights, double mass) const {

    unsigned int fixed_mass = (unsigned int) std::round(mass * WEIGHT_SELECTION_SCALER);
    return (weights.find(fixed_mass) != weights.end());
}

void FragmentGraph::getRandomSampledTransitions(std::set<int> &selected_trans_id, int max_selection) {

    std::vector<int> ids(transitions.size());
    std::iota(ids.begin(), ids.end(), 0);
    std::shuffle(ids.begin(), ids.end(), util_rng);

    max_selection = std::min(max_selection, (int)transitions.size());
    for (int i = 0; i < max_selection; ++i)
        selected_trans_id.insert(ids[i]);
}

void FragmentGraph::clearAllSmiles() {
    for (auto &fragment: fragments) {
        fragment->clearSmiles();
    }
};

//Direct constructor that bipasses the mols altogether and directly sets the nl_smiles
int EvidenceFragmentGraph::addToGraphDirectNoCheck(const EvidenceFragment &fragment, const Transition *transition,
                                                   int parentid) {

    auto id = fragments.size();
    fragments.push_back(EvidenceFragment(fragment, id));
    from_id_tmap.resize(id + 1);
    to_id_tmap.resize(id + 1);
    if (parentid >= 0) addTransition(parentid, id, transition->getNLSmiles());
    return id;
}

void EvidenceFragmentGraph::addTransition(int from_id, int to_id, const std::string *nl_smiles) {
    auto idx = transitions.size();
    transitions.push_back(std::make_shared<Transition>(from_id, to_id, nl_smiles));
    from_id_tmap[from_id].push_back(idx);
    to_id_tmap[to_id].push_back(idx);
}

void EvidenceFragmentGraph::writeFragmentsOnly(std::ostream &out) const {

    for (auto &it : fragments) {
        out << it.getId() << " ";
        out << std::setprecision(10) << it.getMass() << " ";
        out << *(it.getIonSmiles()) << std::endl;
    }
}

void EvidenceFragmentGraph::writeFragmentsOnlyForIds(std::ostream &out, std::set<int> & ids) const {

    for (auto &it : fragments) {
        if (ids.find(it.getId()) != ids.end()) {
            out << it.getId() << " ";
            out << std::setprecision(6) << it.getMass() << " ";
            out << *(it.getIonSmiles()) << std::endl;
        }
    }
}

//Write the FragmentGraph to file (formerly the transition output - without feature details)
void EvidenceFragmentGraph::writeFullGraph(std::ostream &out) const {

    //Fragments
    out << fragments.size() << std::endl;
    writeFragmentsOnly(out);
    out << std::endl;

    //Transitions
    auto itt = transitions.begin();
    for (; itt != transitions.end(); ++itt) {
        out << (*itt)->getFromId() << " ";
        out << (*itt)->getToId() << " ";
        out << *(*itt)->getNLSmiles() << std::endl;
    }
}

bool EvidenceFragmentGraph::fragmentIsRedundant(unsigned int fidx, std::vector<int> &annotated_flags,
                                                std::vector<int> &direct_flags, std::vector<bool> &visited) const {

    if (annotated_flags[fidx])
        return false;

    visited[fidx] = true;

    auto it = from_id_tmap[fidx].begin();
    for (; it != from_id_tmap[fidx].end(); ++it) {
        int cidx = transitions[*it]->getToId();
        if(visited[cidx])
            continue;
        if (!fragmentIsRedundant(cidx, annotated_flags, direct_flags, visited) && !direct_flags[cidx])
            return false;
    }
    return true;
}

void EvidenceFragmentGraph::setFlagsForDirectPaths(std::vector<int> &direct_flags, unsigned int fidx,
                                                   std::vector<int> &annotated_flags) const {

    if (!annotated_flags[fidx])
        return;

    if (direct_flags[fidx])
        return;

    direct_flags[fidx] = 1;
    auto it = from_id_tmap[fidx].begin();
    for (; it != from_id_tmap[fidx].end(); ++it) {
        int cidx = transitions[*it]->getToId();
        setFlagsForDirectPaths(direct_flags, cidx, annotated_flags);
    }
}
