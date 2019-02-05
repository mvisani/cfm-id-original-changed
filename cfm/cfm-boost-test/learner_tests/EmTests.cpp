/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# param_test.cpp
#
# Description: Test code for Param.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

namespace bdata = boost::unit_test::data;
#include "LearnerTestMols.h"
#include "mpi.h"

#include "MolData.h"
#include "EmNNModel.h"
#include "LearnerTestMols.h"

#include <boost/filesystem.hpp>
#include <GraphMol/RDKitBase.h>


void compareSpectra(const Spectrum *orig_spec, const Spectrum *predicted_spec, double abs_mass_tol, double ppm_mass_tol,
                    double intensity_tol) {

    Spectrum::const_iterator ito = orig_spec->begin();
    for (; ito != orig_spec->end(); ++ito) {

        //Find a peak in the predicted spectrum with the same mass
        Spectrum::const_iterator itp = predicted_spec->begin();
        bool found = false;
        for (; itp != predicted_spec->end(); ++itp) {
            double mass_tol = getMassTol(abs_mass_tol, ppm_mass_tol, itp->mass);
            if (fabs(itp->mass - ito->mass) < mass_tol) {
                BOOST_CHECK_SMALL(itp->intensity - ito->intensity,intensity_tol);
                found = true;
                break;
            }
        }
        BOOST_TEST(found);
    }
}

BOOST_AUTO_TEST_SUITE(EMTests)

    BOOST_AUTO_TEST_CASE(EMSelfProdctuctionTest) {
        //set up mpi
        int mp_init_flag;
        MPI_Initialized(&mp_init_flag);
        if(mp_init_flag == 0)
            MPI::Init();

        double intensity_tol = 3.0; //For intensity in range 0-100

        int mpi_rank, mpi_nump;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);

        BOOST_CHECK_EQUAL(mpi_nump, 1);

        //Config
        config_t orig_cfg;
        std::string param_cfg_file = "test_data/example_param_config.txt";
        initConfig(orig_cfg, param_cfg_file);
        orig_cfg.lambda = 0.0000001;
        orig_cfg.use_single_energy_cfm = 1;
        orig_cfg.spectrum_depths[1] = 2;
        orig_cfg.spectrum_depths[2] = 2;
        orig_cfg.include_h_losses = true;

        //Feature Calculator
        std::string feature_cfg_file = "test_data/example_feature_config_withquadratic.txt";
        FeatureCalculator fc(feature_cfg_file);

        //Prepare some simple data
        std::vector<MolData> data;
        std::string id = "TestMol", smiles = "NCCCN";
        std::string spec_file = "test_data/example_spectra.txt";
        data.push_back(MolData(id, smiles, 0, &orig_cfg));
        data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);
        data[0].readInSpectraFromFile(spec_file);

        Param *final_params;
        for (int energy = 0; energy < 3; energy++) {

            config_t cfg;
            initSingleEnergyConfig(cfg, orig_cfg, energy);

            //Run EM
            std::string status_file = "tmp_status_file.log";
            std::string tmp_file = "tmp.log";
            EmModel em(&cfg, &fc, status_file);
            em.trainModel(data, 1, tmp_file, energy);
            std::string param_filename = "tmp_param_output.log";
            em.writeParamsToFile(param_filename);

            if (energy == 0) final_params = new Param(param_filename);
            else {
                Param eparam(param_filename);
                final_params->appendNextEnergyParams(eparam, energy);
            }
        }

        //Predict the output spectra
        data[0].computePredictedSpectra(*final_params);
        data[0].postprocessPredictedSpectra(100.0, 0, 1000);

        //Compare the original and predicted spectra - should be able to overfit
        //very close to the actual values since training on same (and only same) mol
        for (unsigned int energy = 0; energy < data[0].getNumSpectra(); energy++) {
            const Spectrum *orig_spec = data[0].getSpectrum(energy);
            const Spectrum *predicted_spec = data[0].getPredictedSpectrum(energy);
            compareSpectra(orig_spec, predicted_spec, orig_cfg.abs_mass_tol, orig_cfg.ppm_mass_tol, intensity_tol);
        }
    }

    BOOST_AUTO_TEST_CASE(EmNNSelfProdctuctionTest) {
        //set up mpi
        int mp_init_flag;
        MPI_Initialized(&mp_init_flag);
        if(mp_init_flag == 0)
            MPI::Init();

        double intensity_tol = 3.0; //For intensity in range 0-100

        int mpi_rank, mpi_nump;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);

        BOOST_CHECK_EQUAL(mpi_nump, 1);

        //Config
        config_t orig_cfg;
        std::string param_cfg_file = "test_data/example_nnparam_config.txt";
        initConfig(orig_cfg, param_cfg_file);
        orig_cfg.lambda = 0.0000001;
        orig_cfg.use_single_energy_cfm = 1;
        //Feature Calculator
        std::string feature_cfg_file = "test_data/example_feature_config_withquadratic.txt";
        FeatureCalculator fc(feature_cfg_file);

        //Prepare some simple data
        std::vector<MolData> data;
        std::string id = "TestMol", smiles = "NCCCN";
        std::string spec_file = "test_data/example_spectra.txt";
        data.push_back(MolData(id, smiles, 0, &orig_cfg));
        data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);
        data[0].readInSpectraFromFile(spec_file);

        NNParam *final_params;
        for (int energy = 0; energy < 3; energy++) {

            config_t cfg;
            initSingleEnergyConfig(cfg, orig_cfg, energy);
            std::string param_filename = "tmp_param_output.log";

            //Run EM (multiple times, and take best Q)
            double best_Q = -1000000.0;
            const int trials_max = 1;
            std::vector<double> Qs;
            for (int trial = 0; trial < trials_max; trial++) {
                std::string status_file = "tmp_status_file.log";
                std::string tmp_file = "tmp.log";
                EmNNModel em(&cfg, &fc, status_file);
                double Q = em.trainModel(data, 1, tmp_file, energy);
                Qs.push_back(Q);

                if (Q > best_Q) {
                    em.writeParamsToFile(param_filename);
                    best_Q = Q;
                }
            }
            std::vector<double>::iterator it = Qs.begin();
            for (; it != Qs.end(); ++it) std::cout << *it << " ";
            std::cout << " Best=" << best_Q << std::endl;

            if (energy == 0)
                final_params = new NNParam(param_filename);
            else {
                NNParam eparam(param_filename);
                final_params->appendNextEnergyParams(eparam, energy);
            }
        }

        //Predict the output spectra
        data[0].computePredictedSpectra(*final_params);
        data[0].postprocessPredictedSpectra(100.0, 0);

        //Compare the original and predicted spectra - should be able to overfit
        //very close to the actual values since training on same (and only same) mol
        for (unsigned int energy = 0; energy < data[0].getNumSpectra(); energy++) {
            const Spectrum *orig_spec = data[0].getSpectrum(energy);
            const Spectrum *predicted_spec = data[0].getPredictedSpectrum(energy);
            std::cout << energy << " ";
            compareSpectra(orig_spec, predicted_spec, orig_cfg.abs_mass_tol, orig_cfg.ppm_mass_tol, intensity_tol);
        }
    }

    BOOST_AUTO_TEST_CASE(EMSelfProdctuctionWithIsotopeTest) {
        //set up mpi
        int mp_init_flag;
        MPI_Initialized(&mp_init_flag);
        if(mp_init_flag == 0)
            MPI::Init();

        double intensity_tol = 2.0; //For intensity in range 0-100

        int mpi_rank, mpi_nump;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);

        BOOST_CHECK_EQUAL(mpi_nump, 1);

        //Config
        config_t orig_cfg;
        std::string param_cfg_file = "test_data/example_param_config.txt";
        initConfig(orig_cfg, param_cfg_file);
        orig_cfg.lambda = 0.0000001;
        orig_cfg.use_single_energy_cfm = 1;
        orig_cfg.spectrum_depths.resize(1);
        orig_cfg.spectrum_weights.resize(1);
        orig_cfg.include_isotopes = 1;
        orig_cfg.ga_converge_thresh = 0.0001;
        orig_cfg.disable_training_metrics = true;

        //Feature Calculator
        std::string feature_cfg_file = "test_data/example_feature_config_withquadratic.txt";
        FeatureCalculator fc(feature_cfg_file);

        ///Prepare some simple data
        std::vector<MolData> data;
        data.push_back(IsotopeTestMol(&orig_cfg));
        data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, true);

        Param *final_params;

        config_t cfg;
        initSingleEnergyConfig(cfg, orig_cfg, 0);

        //Run EM
        std::string status_file = "tmp_status_file.log";
        std::string tmp_file = "tmp.log";
        EmModel em(&cfg, &fc, status_file);
        em.trainModel(data, 1, tmp_file, 0);
        std::string param_filename = "tmp_param_output.log";
        em.writeParamsToFile(param_filename);
        final_params = new Param(param_filename);

        //Predict the output spectra
        data[0].computePredictedSpectra(*final_params);
        data[0].postprocessPredictedSpectra(100.0, 0, 1000);

        const Spectrum *orig_spec = data[0].getSpectrum(0);
        const Spectrum *predicted_spec = data[0].getPredictedSpectrum(0);
        compareSpectra(orig_spec, predicted_spec, orig_cfg.abs_mass_tol, orig_cfg.ppm_mass_tol, intensity_tol);

    }

    BOOST_AUTO_TEST_CASE(EMTestMiniBatchSelection) {

        //set up mpi
        int mp_init_flag;
        MPI_Initialized(&mp_init_flag);
        if(mp_init_flag == 0)
            MPI::Init();

        config_t cfg;
        initDefaultConfig(cfg);
        cfg.ga_minibatch_nth_size = 10; //Take 1 in 10
        initDerivedConfig(cfg);

        std::string status_filename = "dummy_status.log";
        std::vector<std::string> fnames;
        fnames.push_back("HydrogenRemoval");
        FeatureCalculator fc(fnames);

        std::vector<int> flags1(1000);    //Select 100
        std::vector<int> flags2(1000);

        EmModel em(&cfg, &fc, status_filename);
        em.setMiniBatchFlags(flags1, cfg.ga_minibatch_nth_size);
        em.setMiniBatchFlags(flags2, cfg.ga_minibatch_nth_size);

        //Check that the random selections select the right number of molecules,
        //and that the selected molecules are different in the two runs
        int num_on1 = std::count(flags1.begin(), flags1.end(), 1);
        int num_on2 = std::count(flags2.begin(), flags2.end(), 1);
        BOOST_CHECK_EQUAL(num_on1, 100);
        BOOST_CHECK_EQUAL(num_on1, 100);

        std::vector<int>::iterator it1 = flags1.begin(), it2 = flags2.begin();
        int count_same = 0;
        for (; it1 != flags1.end(); ++it1, ++it2) {
            count_same += (*it1) * (*it2);
        }
        BOOST_TEST(count_same != num_on1 );
        BOOST_TEST(count_same != num_on2 );
    }
BOOST_AUTO_TEST_SUITE_END()
