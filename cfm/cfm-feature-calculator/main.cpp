/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Learn the parameters of a model for the mass spec
#				fragmentation process using EM, then use it to
predict #				the mass spectra.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "mpi.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "Config.h"
#include "EmModel.h"

void parseInputFile(std::vector<std::string> &data, std::string &input_filename,
                    int mpi_rank, int mpi_nump);

int main(int argc, char *argv[]) {
    int mpi_rank, mpi_nump;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);
    MPI_Status status;

    if (argc != 5) {
        std::cout << std::endl << std::endl;
        std::cout
                << std::endl
                << "Usage: cfm-feature-calculator <input_filename> <feature_filename> "
                << std::endl
                << std::endl
                << std::endl;
        std::cout << std::endl
                  << "input_filename:" << std::endl
                  << "Text file with number of mols on first line, then "
                  << std::endl
                  << "id smiles_or_inchi cross_validation_group" << std::endl
                  << "on each line after that." << std::endl;
        std::cout
                << std::endl
                << "feature_filename:" << std::endl
                << "Text file with list of feature names to include, line separated:"
                << std::endl
                << "BreakAtomPair" << std::endl
                << "IonRootPairs...etc" << std::endl;
        std::cout << std::endl << "save_filename: csv file to save ouput" << std::endl;
        exit(1);
    }

    std::string input_filename = argv[1]; // List (one per line): id, smiles_or_inchi, group
    std::string feature_filename = argv[2]; // List of features, line-spaced
    std::string config_filename = argv[3];  // Parameter configuration
    std::string save_filename = argv[4];    // MSP file or Directory containing
    // the peak files for each molecule
    // (in format <id>.txt)

    std::string status_filename("status.log"); // Status file to write to

    MPI_Barrier(MPI_COMM_WORLD);

    // Init feature calculator
    if (mpi_rank == MASTER)
        std::cout << "Initialising Feature Calculator..";
    FeatureCalculator fc(feature_filename);
    if (mpi_rank == MASTER)
    {
        for(auto featureName : fc.getFeatureNames())
        {
            std::cout << featureName << std::endl;
        }
        std::cout << "Total number of features: " << fc.getNumFeatures()  << std::endl;
        std::cout << "Done" << std::endl;
    }

    // Init config file
    if (mpi_rank == MASTER)
        std::cout << "Initialising Parameter Configuration..";
    config_t cfg;
    initConfig(cfg, config_filename, mpi_rank == MASTER);
    if (mpi_rank == MASTER)
        std::cout << "Done" << std::endl;

    // Read Mols from file
    if (mpi_rank == MASTER)
        std::cout << "Parsing input file...";
    /*std::vector<MolData> data;
    parseInputFile(data, input_filename, mpi_rank, mpi_nump, &cfg);*/
    std::vector<std::string> data;
    parseInputFile(data, input_filename, mpi_rank, mpi_nump);

    if (mpi_rank == MASTER)
        std::cout << "Done" << std::endl;
    // MPI_Barrier(MPI_COMM_WORLD);

    // Fragment Graph Computation (or load from file)
    time_t before_fg, after_fg;
    before_fg = time(nullptr);
    if (mpi_rank == MASTER)
        std::cout << "Computing fragmentation graphs and features..";

    MPI_File output_file;
    MPI_File_open(MPI_COMM_WORLD, save_filename.c_str(),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &output_file);

    int success_count = 0, except_count = 0;

    for (auto molStr : data) {

        std::stringstream ss(molStr);
        std::string id;
        std::string smiles_or_inchi;
        int group;
        ss >> id >> smiles_or_inchi >> group;
        MolData *mol = new MolData(id, smiles_or_inchi, group, &cfg);

        try {
            time_t before, after;
            before = time(nullptr);
            mol->ComputeFragmentGraphAndReplaceMolsWithFVs(&fc);
            std::ofstream eout;
            eout.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
            after = time(nullptr);
            eout << mol->getId() << "Done. Time Elaspsed = " << (after - before)
                 << " Seconds";
            eout << " :Num Frag = " << mol->GetFragmentGraph()->getNumFragments();
            eout << " :Num Trans = " << mol->GetFragmentGraph()->getNumTransitions()
                 << std::endl;
            eout.close();

            //std::string molFeatureCsvStr = mol->getFVsAsCSVString();
            std::string molFeatureCsvStr = mol->GetFVsAsSparseCSVString();
            MPI_File_write_shared(output_file, molFeatureCsvStr.c_str(), molFeatureCsvStr.size(), MPI_CHAR, &status);

            success_count++;
            delete mol;
        }
        catch (std::exception &e) {
            std::ofstream eout;
            eout.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
            eout << "Exception occurred computing fragment graph for " << mol->getId()
                 << std::endl;
            eout << mol->getSmilesOrInchi() << std::endl;
            eout << e.what() << std::endl << std::endl;
            except_count++;
            eout << except_count << " exceptions, from "
                 << except_count + success_count << " total" << std::endl;
            eout.close();
            delete mol;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    after_fg = time(nullptr);
    if (mpi_rank == MASTER)
        std::cout << "Done" << std::endl;
    std::cout << mpi_rank << ": " << success_count << " successfully computed. "
              << except_count << " exceptions." << std::endl;
    if (mpi_rank == MASTER) {
        std::cout << "Total Fragmentation Graph Computation Time Elaspsed = "
                  << (after_fg - before_fg) << " Seconds" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);    //Wait for all threads
    MPI_File_close(&output_file);
    MPI_Finalize();
    return (0);

}

void parseInputFile(std::vector<std::string> &data, std::string &input_filename,
                    int mpi_rank, int mpi_nump) {

    std::string line, smiles_or_inchi, id;
    std::ifstream ifs(input_filename.c_str(), std::ifstream::in);
    int num_mols = 0;

    // Get the first line - the number of input molecules
    if (ifs.good()) {
        getline(ifs, line);
        num_mols = atoi(line.c_str());
        if (mpi_rank == MASTER)
            std::cout << "Reading " << num_mols << " mols" << std::endl;
    } else {
        std::cout << "Could not open input file " << input_filename << std::endl;
    }

    // Now get all the molecules
    int i = 0;
    int idx = 0;
    while (ifs.good() && i < num_mols) {
        i++;
        getline(ifs, line);
        if (line.size() < 3)
            continue;
        // split data into each processor
        // each processor should have even amount of mols
        idx++;
        if ((i % mpi_nump) == mpi_rank)
            data.push_back(line);
    }
    std::cout << mpi_rank << ":Data Size: " << data.size() << std::endl;
}
