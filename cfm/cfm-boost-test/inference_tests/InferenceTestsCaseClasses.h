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

#include "Config.h"
#include "FragmentTreeNode.h"
#include "MolData.h"

#include <boost/filesystem.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

#pragma onceonce

// Fixture requires more code to do the same thing
// since we only have one test case, let us just use this
class InferTestMol : public MolData {
public:
    InferTestMol(config_t *cfg) : MolData("Infer Test Mol", "", cfg) {

        //Create a molecule based on what was previously in test_bn_transition_ipfp.txt
        fg = new FragmentGraph();
        FeatureHelper fh;
        romol_ptr_t basic_nl = createMolPtr("C");
        std::vector<int> null_eloc;

        fg->addToGraph(FragmentTreeNode(createMolPtr("O=C(O)C[NH2+]C(=O)C(NC(=O)CN)Cc1ccccc1"), basic_nl, -1, -1, &fh,
                                        null_eloc), -1); //id = 0
        fg->addToGraph(
                FragmentTreeNode(createMolPtr("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1"), basic_nl, -1, -1, &fh, null_eloc),
                0); // id = 1, 0 -> 1
        fg->addToGraph(FragmentTreeNode(createMolPtr("N=CC(=O)[NH2+]CCc1ccccc1"), basic_nl, -1, -1, &fh, null_eloc),
                       0); //id = 2, 0 -> 2
        fg->addToGraph(FragmentTreeNode(createMolPtr("[NH2+]=CCc1ccccc1"), basic_nl, -1, -1, &fh, null_eloc),
                       2); // id = 3, 2 -> 3

        //Set thetas/transition probs to match matlab reference2
        thetas.resize(3);
        for (int energy = 0; energy < 3; energy++) {
            thetas[energy].resize(fg->getNumTransitions());
            thetas[energy][0] = 1.386294361119891;    //0->1
            thetas[energy][1] = 2.456735772821304;    //0->2
            thetas[energy][2] = 0.0;    //2->3
        }
        computeLogTransitionProbabilities();

        //Set the spectra
        spectra.resize(3);
        spectra[0].push_back(Peak(120.082039, 35.914500));
        spectra[0].push_back(Peak(177.103613, 50.000000));
        spectra[0].push_back(Peak(280.133689, 50.146650));
        spectra[1].push_back(Peak(120.083820, 100.00000));
        spectra[1].push_back(Peak(177.106284, 33.000500));
        spectra[2].push_back(Peak(120.081802, 100.00000));
        for (int i = 0; i <= 2; i++)
            spectra[i].normalizeAndSort();

    }
};

class InferComplexTestMol : public MolData {
public:
    InferComplexTestMol(config_t *cfg) : MolData("NIST2011_1201", "CC(C)C(N=C(O)CN)C(=O)O", cfg) {

        //Compulte a fragmentation graph
        std::vector<std::string> fnames;
        fnames.push_back("BreakAtomPair");
        FeatureCalculator fc(fnames);
        computeFragmentGraphAndReplaceMolsWithFVs(&fc);

        //Randomly set the theta values
        thetas.resize(1);
        thetas[0].resize(fg->getNumTransitions());
        for (int i = 0; i < fg->getNumTransitions(); i++)
            thetas[0][i] = (double(std::rand()) / double(RAND_MAX) - 0.5) * 3;
        computeLogTransitionProbabilities();

        //Load a very simple spectrum
        spectra.resize(1);
        spectra[0].push_back(Peak(114.0, 100.0));
        spectra[0].push_back(Peak(174.0, 100.0));
        spectra[0].normalizeAndSort();

    }
};

class SpectrumTestMol : public MolData {
public:
    SpectrumTestMol(config_t *cfg) : MolData("Spectrum Test Mol", "", cfg) {

        //Create a molecule based on what was previously in test_bn_transition_ipfp.txt
        fg = new FragmentGraph();
        FeatureHelper fh;
        romol_ptr_t basic_nl = createMolPtr("C");
        //Root node with Mass1 (0), two fragments with Mass2 (1,2), three fragments with Mass3 (3,4,5)
        std::vector<int> null_eloc;
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCCCCCO"), basic_nl, -1, -1, &fh, null_eloc), -1);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CC(C)C"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCCC"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCOC"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CC(C)O"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCCO"), basic_nl, -1, -1, &fh, null_eloc), 0);

        //Set the spectra
        spectra.resize(1);
        spectra[0].push_back(Peak(102.10446506800000, 10.0));
        spectra[0].push_back(Peak(58.078250319999995, 30.0));
        spectra[0].push_back(Peak(60.057514875999999, 60.0));
        spectra[0].normalizeAndSort();

    }
};

class SpectrumTestMolNoisePeak : public MolData {
public:
    SpectrumTestMolNoisePeak(config_t *cfg) : MolData("Spectrum Test Mol with Noise Peak", "", cfg) {

        //Create a molecule  as above, but add a peak quite far from any fragment mass (noise peak)
        fg = new FragmentGraph();
        FeatureHelper fh;
        romol_ptr_t basic_nl = createMolPtr("C");
        std::vector<int> null_eloc;
        //Root node with Mass1 (0), two fragments with Mass2 (1,2), three fragments with Mass3 (3,4,5)
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCCCCCO"), basic_nl, -1, -1, &fh, null_eloc), -1);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CC(C)C"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCCC"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCOC"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CC(C)O"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCCO"), basic_nl, -1, -1, &fh, null_eloc), 0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("CCOO"), basic_nl, -1, -1, &fh, null_eloc), 0);

        //Set the spectra
        spectra.resize(1);
        spectra[0].push_back(Peak(102.10446506800000, 10.0));
        spectra[0].push_back(Peak(58.078250319999995, 30.0));
        spectra[0].push_back(Peak(60.057514875999999, 60.0));
        spectra[0].push_back(Peak(70.0, 50.0));
        spectra[0].normalizeAndSort();

    }
};

class IsotopeSpectrumTestMolNoisePeak : public MolData {
public:
    IsotopeSpectrumTestMolNoisePeak(config_t *cfg) : MolData("Isotope Spectrum Test Mol with Noise Peak", "", cfg) {

        //Create a molecule as for isotope test above, but add a peak quite far from any fragment mass (noise peak)
        fg = new FragmentGraph(cfg);
        FeatureHelper fh;
        romol_ptr_t basic_nl = createMolPtr("C");
        //Two molecules, (almost) same mass but different isotope patterns.
        std::vector<int> null_eloc;
        fg->addToGraph(
                FragmentTreeNode(createMolPtr("[CH2-]C(C)NC1=CN=CC(=C1C#N)Cl"), basic_nl, -1, -1, &fh, null_eloc), -1);
        fg->addToGraph(
                FragmentTreeNode(createMolPtr("[CH-](C1=NNN=N1)NC2=NC(=O)C(=O)N2"), basic_nl, -1, -1, &fh, null_eloc),
                0);
        fg->addToGraph(FragmentTreeNode(createMolPtr("[CH-](C1=NNN=N1)NC"), basic_nl, -1, -1, &fh, null_eloc), 0);

        //Set the spectra
        spectra.resize(1);
        spectra[0].push_back(Peak(194.0490426799, 100.0000000000));
        spectra[0].push_back(Peak(195.0518167914, 11.3124853142));
        spectra[0].push_back(Peak(196.0462415099, 32.9809130321));
        spectra[0].push_back(Peak(197.0489073693, 3.6831305376));
        spectra[0].push_back(Peak(50.0, 100.0));
        spectra[0].normalizeAndSort();

    }
};

class IsotopeSpectrumTestMol : public MolData {
public:
    IsotopeSpectrumTestMol(config_t *cfg) : MolData("Isotope Spectrum Test Mol", "", cfg) {

        //Create a molecule based on what was previously in test_bn_transition_ipfp.txt
        fg = new FragmentGraph(cfg);
        FeatureHelper fh;
        romol_ptr_t basic_nl = createMolPtr("C");
        std::vector<int> null_eloc;
        //Two molecules, (almost) same mass but different isotope patterns.
        fg->addToGraph(
                FragmentTreeNode(createMolPtr("[CH2-]C(C)NC1=CN=CC(=C1C#N)Cl"), basic_nl, -1, -1, &fh, null_eloc), -1);
        fg->addToGraph(
                FragmentTreeNode(createMolPtr("[CH-](C1=NNN=N1)NC2=NC(=O)C(=O)N2"), basic_nl, -1, -1, &fh, null_eloc),
                0);

        //Set the spectra
        spectra.resize(1);
        spectra[0].push_back(Peak(194.0490426799, 100.0000000000));
        spectra[0].push_back(Peak(195.0518167914, 11.3124853142));
        spectra[0].push_back(Peak(196.0462415099, 32.9809130321));
        spectra[0].push_back(Peak(197.0489073693, 3.6831305376));
        spectra[0].normalizeAndSort();

    }
};

