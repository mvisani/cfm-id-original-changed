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

#include <fstream>
#include <iostream>
#include <sstream>

#include "Config.h"
#include "EM.h"
#include "EM_NN.h"
#include "Feature.h"
#include "FragmentGraphGenerator.h"
#include "MolData.h"

void parseInputFile(std::vector<MolData> &data, std::string &input_filename,
                    int mpi_rank, int mpi_nump, config_t *cfg);

int main(int argc, char *argv[]) {
  int mpi_rank, mpi_nump;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);

  if (argc != 4) {
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
    std::cout << std::endl
              << "output_filename:" << std::endl
              << "CSV file used as output" << std::endl std::cout << std::endl;
    exit(1);
  }

  std::string input_filename =
      argv[1]; // List (one per line): id, smiles_or_inchi, group
  std::string feature_filename = argv[2]; // List of features, line-spaced
  std::string config_filename = argv[3];  // Parameter configuration
  std::string save_filename = argv[4];    // MSP file or Directory containing
                                          // the peak files for each molecule
                                          // (in format <id>.txt)

  std::string status_filename("status.log"); // Status file to write to

  if (mpi_rank == MASTER) {
    // Create the tmp_data directory if it doesn't exist
    if (!boost::filesystem::exists("tmp_data"))
      boost::filesystem::create_directory("tmp_data");
    if (!boost::filesystem::exists("tmp_data/enumerated_output"))
      boost::filesystem::create_directory("tmp_data/enumerated_output");
    if (!boost::filesystem::exists("tmp_data/predicted_output"))
      boost::filesystem::create_directory("tmp_data/predicted_output");
    if (!boost::filesystem::exists("tmp_data/fv_fragment_graphs"))
      boost::filesystem::create_directory("tmp_data/fv_fragment_graphs");
    // Delete the status file if it already exists
    if (boost::filesystem::exists(status_filename))
      boost::filesystem::remove_all(status_filename);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (mpi_rank == MASTER)
    std::cout << "Initialising Feature Calculator..";
  FeatureCalculator fc(feature_filename);
  if (mpi_rank == MASTER)
    std::cout << "Done" << std::endl;

  if (mpi_rank == MASTER)
    std::cout << "Initialising Parameter Configuration..";
  config_t cfg;
  initConfig(cfg, config_filename, mpi_rank == MASTER);
  if (mpi_rank == MASTER)
    std::cout << "Done" << std::endl;

  if (mpi_rank == MASTER)
    std::cout << "Parsing input file...";
  std::vector<MolData> data;
  parseInputFile(data, input_filename, mpi_rank, mpi_nump, &cfg);
  if (mpi_rank == MASTER)
    std::cout << "Done" << std::endl;

  // Configure fragment graph state files
  std::string fv_filename_out = "tmp_data/fv_fragment_graphs/P" +
                                boost::lexical_cast<std::string>(mpi_rank) +
                                "_graphs.fg";
  std::string fv_filename_in = "tmp_data/fv_fragment_graphs/P" +
                               boost::lexical_cast<std::string>(mpi_rank) +
                               "_graphs.fg_tmp";
  if (min_group != 0)
    fv_filename_in = fv_filename_out;
  if (min_group == 0 && boost::filesystem::exists(fv_filename_in))
    boost::filesystem::remove(fv_filename_in);
  if (min_group == 0 && boost::filesystem::exists(fv_filename_out))
    boost::filesystem::copy_file(fv_filename_out, fv_filename_in);
  std::ifstream *fv_ifs;
  std::string next_id = "";
  if (boost::filesystem::exists(fv_filename_in)) {
    fv_ifs = new std::ifstream(fv_filename_in.c_str(),
                               std::ifstream::in | std::ios::binary);
    if (!(*fv_ifs))
      std::cout << "Could not open file " << fv_filename_in << std::endl;
    else {
      unsigned int id_size;
      fv_ifs->read(reinterpret_cast<char *>(&id_size), sizeof(id_size));
      next_id.resize(id_size);
      fv_ifs->read(&next_id[0], id_size);
    }
  }
  std::ofstream fv_out;
  if (min_group == 0)
    fv_out.open(fv_filename_out.c_str(), std::ios::out | std::ios::binary);

  // Fragment Graph Computation (or load from file)
  time_t before_fg, after_fg;
  before_fg = time(NULL);
  if (mpi_rank == MASTER)
    std::cout << "Computing fragmentation graphs and features..";
  std::vector<MolData>::iterator mit = data.begin();
  int success_count = 0, except_count = 0;
  for (; mit != data.end(); ++mit) {
    try {
      // If we're not training, only load the ones we'll be testing
      if ((mit->getGroup() >= min_group && mit->getGroup() <= max_group) ||
          !no_train) {
        if (mit->getId() == next_id) {
          mit->readInFVFragmentGraphFromStream(*fv_ifs);
          unsigned int id_size;
          fv_ifs->read(reinterpret_cast<char *>(&id_size), sizeof(id_size));
          if (!fv_ifs->eof()) {
            next_id.resize(id_size);
            fv_ifs->read(&next_id[0], id_size);
          } else
            next_id = "NULL_ID";

        } else {
          time_t before, after;
          before = time(NULL);
          mit->computeFragmentGraphAndReplaceMolsWithFVs(&fc);
          std::ofstream eout;
          eout.open(status_filename.c_str(),
                    std::fstream::out | std::fstream::app);
          after = time(NULL);
          eout << mit->getId() << "Done. Time Elaspsed = " << (after - before)
               << " Seconds";
          eout << " :Num Frag = " << mit->getFragmentGraph()->getNumFragments();
          eout << " :Num Trans = "
               << mit->getFragmentGraph()->getNumTransitions() << std::endl;
          eout.close();
        }
        unsigned int id_size = mit->getId().size();
        if (min_group == 0) {
          fv_out.write(reinterpret_cast<const char *>(&id_size),
                       sizeof(id_size));
          fv_out.write(&(mit->getId()[0]), id_size);
          mit->writeFVFragmentGraphToStream(fv_out); // We always write it, in
                                                     // case we haven't already
                                                     // computed all of them
        }
        success_count++;
      }
    } catch (std::exception e) {
      std::ofstream eout;
      eout.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
      eout << "Exception occurred computing fragment graph for " << mit->getId()
           << std::endl;
      eout << mit->getSmilesOrInchi() << std::endl;
      eout << e.what() << std::endl << std::endl;
      except_count++;
      eout << except_count << " exceptions, from "
           << except_count + success_count << " total" << std::endl;
      eout.close();
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  after_fg = time(NULL);
  if (mpi_rank == MASTER)
    std::cout << "Done" << std::endl;
  std::cout << mpi_rank << ": " << success_count << " successfully computed. "
            << except_count << " exceptions." << std::endl;
  if (mpi_rank == MASTER)
    std::cout << "Total Fragmentation Graph Computation Time Elaspsed = "
              << (after_fg - before_fg) << " Seconds";
  if (min_group == 0)
    fv_out.close();
  if (boost::filesystem::exists(fv_filename_in)) {
    if (*fv_ifs)
      fv_ifs->close();
    delete fv_ifs;
    if (min_group == 0)
      boost::filesystem::remove(fv_filename_in);
  }
  MPI_Finalize();

  // Fragment Graph Computation (or load from file)
  time_t before_fg, after_fg;
  before_fg = time(NULL);
  if (mpi_rank == MASTER)
    std::cout << "Computing fragmentation graphs and features..";
  std::vector<MolData>::iterator mit = data.begin();
  int success_count = 0, except_count = 0;
  for (; mit != data.end(); ++mit) {
    try {
      // If we're not training, only load the ones we'll be testing
      if ((mit->getGroup() >= min_group && mit->getGroup() <= max_group) ||
          !no_train) {
        if (mit->getId() == next_id) {
          mit->readInFVFragmentGraphFromStream(*fv_ifs);
          unsigned int id_size;
          fv_ifs->read(reinterpret_cast<char *>(&id_size), sizeof(id_size));
          if (!fv_ifs->eof()) {
            next_id.resize(id_size);
            fv_ifs->read(&next_id[0], id_size);
          } else
            next_id = "NULL_ID";

        } else {
          time_t before, after;
          before = time(NULL);
          mit->computeFragmentGraphAndReplaceMolsWithFVs(&fc);
          std::ofstream eout;
          eout.open(status_filename.c_str(),
                    std::fstream::out | std::fstream::app);
          after = time(NULL);
          eout << mit->getId() << "Done. Time Elaspsed = " << (after - before)
               << " Seconds";
          eout << " :Num Frag = " << mit->getFragmentGraph()->getNumFragments();
          eout << " :Num Trans = "
               << mit->getFragmentGraph()->getNumTransitions() << std::endl;
          eout.close();
        }
        unsigned int id_size = mit->getId().size();
        if (min_group == 0) {
          fv_out.write(reinterpret_cast<const char *>(&id_size),
                       sizeof(id_size));
          fv_out.write(&(mit->getId()[0]), id_size);
          mit->writeFVFragmentGraphToStream(fv_out); // We always write it, in
                                                     // case we haven't already
                                                     // computed all of them
        }
        success_count++;
      }
    } catch (std::exception e) {
      std::ofstream eout;
      eout.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
      eout << "Exception occurred computing fragment graph for " << mit->getId()
           << std::endl;
      eout << mit->getSmilesOrInchi() << std::endl;
      eout << e.what() << std::endl << std::endl;
      except_count++;
      eout << except_count << " exceptions, from "
           << except_count + success_count << " total" << std::endl;
      eout.close();
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  after_fg = time(NULL);
  if (mpi_rank == MASTER)
    std::cout << "Done" << std::endl;
  std::cout << mpi_rank << ": " << success_count << " successfully computed. "
            << except_count << " exceptions." << std::endl;
  if (mpi_rank == MASTER)
    std::cout << "Total Fragmentation Graph Computation Time Elaspsed = "
              << (after_fg - before_fg) << " Seconds";
  if (min_group == 0)
    fv_out.close();
  if (boost::filesystem::exists(fv_filename_in)) {
    if (*fv_ifs)
      fv_ifs->close();
    delete fv_ifs;
    if (min_group == 0)
      boost::filesystem::remove(fv_filename_in);
  }

  return (0);
}

void parseInputFile(std::vector<MolData> &data, std::string &input_filename,
                    int mpi_rank, int mpi_nump, config_t *cfg) {

  std::string line, smiles_or_inchi, id;
  std::ifstream ifs(input_filename.c_str(), std::ifstream::in);
  int group, num_mols = 0;

  // Get the first line - the number of input molecules
  if (ifs.good()) {
    getline(ifs, line);
    num_mols = atoi(line.c_str());
  } else {
    std::cout << "Could not open input file " << input_filename << std::endl;
  }

  // Now get all the molecules
  int i = 0;
  while (ifs.good() && i < num_mols) {
    i++;

    getline(ifs, line);
    if (line.size() < 3)
      continue;

    std::stringstream ss(line);
    ss >> id >> smiles_or_inchi >> group;

    // Split the data between processors. Only load in data for this
    // processor
    if ((i % mpi_nump) == mpi_rank)
      data.push_back(MolData(id, smiles_or_inchi, group, cfg));
  }
}
