/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Compute spectrum prediction summary statistics.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Comparators.h"

#include <GraphMol/SanitException.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

void reportMeanStd(std::ostream &out, std::vector<double> &scores);

void parseInputFile(std::vector<MolData> &data, std::string &input_filename);

int main(int argc, char *argv[]) {
    bool to_stdout = true;
    std::string output_filename, input_filename;
    std::string measured_spec_dir, predicted_spec_dir;
    double ppm_mass_tol = 10.0, abs_mass_tol = 0.01;
    int num_spectra = 0;

    if (argc < 5 || argc > 14) {
        std::cout << std::endl << std::endl;
        std::cout << std::endl
                  << "Usage: compute-stats <input_filename> "
                     "<measured_spec_dir> <predicted_spec_dir> "
                     "<num_spectra_per_mol> <ppm_mass_tol> <abs_mass_tol> "
                     "<output_filename> <cumulative_intensity_thresh> "
                     "<apply_cutoffs> <clean_target_spectra> "
                     "<quantise_spectra_dec_pl> <num_groups>"
                  << std::endl
                  << std::endl
                  << std::endl;
        std::cout << std::endl
                  << "input_filename:" << std::endl
                  << "Input file containing list of ids to compute stats for (same "
                     "input format as to cfm-train)."
                  << std::endl;
        std::cout << std::endl
                  << "measured_spec_dir:" << std::endl
                  << "Directory containing measured spectra or msp file."
                  << std::endl;
        std::cout << std::endl
                  << "predicted_spec_dir:" << std::endl
                  << "Directory containing predicted spectra or msp file (use P#G$ "
                     "for general msp filename)"
                  << std::endl;
        std::cout << std::endl
                  << "num_spectra_per_mol:" << std::endl
                  << "Number of spectra expected to be found for each molecule"
                  << std::endl;
        std::cout << std::endl
                  << "ppm_mass_tol (opt):" << std::endl
                  << "The mass tolerance in ppm to use when matching peaks within "
                     "the dot product comparison - will use higher resulting "
                     "tolerance of ppm and abs (if not given defaults to 10ppm)"
                  << std::endl;
        std::cout << std::endl
                  << "abs_mass_tol (opt):" << std::endl
                  << "The mass tolerance in abs Da to use when matching peaks "
                     "within the dot product comparison - will use higher "
                     "resulting tolerance of ppm and abs ( if not given defaults "
                     "to 0.01Da)"
                  << std::endl;
        std::cout << std::endl
                  << "output_filename (opt):" << std::endl
                  << "The filename of the output file to write to (if not given, "
                     "prints to stdout)"
                  << std::endl;
        std::cout << std::endl
                  << "cumulative_intensity_thresh (opt):" << std::endl
                  << "Ordering peaks from maximum to minimum intensities, take "
                     "peaks until the sum total of their intensities as a "
                     "percentage of the total is above this value (default = 100 "
                     "(off))"
                  << std::endl;
        std::cout << std::endl
                  << "minimum relative intensity (opt):" << std::endl
                  << "Filter out any peaks below this relative intensity (default = 0, range 0-100)"
                  << std::endl;
        std::cout << std::endl
                  << "minimum peak cutoffs (opt) (default = 1))" << std::endl
                  << std::endl;
        std::cout << std::endl
                  << "maximum peak cutoffs (opt) (default = 1000000):"
                  << std::endl;
        std::cout << std::endl
                  << "clean_target_spectra (opt):" << std::endl
                  << "Whether to clean the target spectra before computing stats "
                     "for them (default = 0(off))"
                  << std::endl;
        std::cout << std::endl
                  << "quantise_spectra_dec_pl (opt):" << std::endl
                  << "Quantise the measured and predicted spectra to masses of "
                     "this number of decimal places (default -1(off))"
                  << std::endl;
        std::cout << std::endl
                  << "group_to_compute (opt):" << std::endl
                  << "Restrict computation to this group only (Default -1(ignore))"
                  << std::endl;
        exit(1);
    }

    input_filename = argv[1];
    int group_to_compute = -1;
    if (argc >= 14) {
        try {
            group_to_compute = boost::lexical_cast<bool>(argv[12]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid group_to_compute (Expecting numerical): "
                      << argv[13] << std::endl;
            exit(1);
        }
    }

    // Spectrum Inputs
    MspReader *measured_msp, *predicted_msp;
    measured_spec_dir = argv[2];
    bool measured_is_msp = false;
    std::string measured_pre_id = "";
    if (measured_spec_dir.size() >= 4 &&
        measured_spec_dir.substr(measured_spec_dir.size() - 4, 4) == ".msp") {
        measured_is_msp = true;
        measured_msp =
                new MspReader(measured_spec_dir.c_str(), measured_pre_id.c_str());
    }
    predicted_spec_dir = argv[3];
    bool predicted_is_msp = false;
    std::string predicted_pre_id = "";
    if (predicted_spec_dir.size() >= 4 &&
        predicted_spec_dir.substr(predicted_spec_dir.size() - 4, 4) == ".msp") {
        predicted_is_msp = true;
        if (group_to_compute != -1) {
            std::string group_msp_file = boost::algorithm::replace_all_copy(
                    predicted_spec_dir, "$",
                    boost::lexical_cast<std::string>(group_to_compute));
            predicted_msp = new MspReader(group_msp_file.c_str(), "");
        } else
            predicted_msp =
                    new MspReader(predicted_spec_dir.c_str(), predicted_pre_id.c_str());
    }

    // Other inputs
    if (argc >= 5) {
        try {
            num_spectra = boost::lexical_cast<int>(argv[4]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid num_spectra (Expecting numerical): " << argv[4]
                      << std::endl;
            exit(1);
        }
    }
    if (argc >= 6) {
        try {
            ppm_mass_tol = boost::lexical_cast<double>(argv[5]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid ppm_mass_tol (Expecting numerical): " << argv[5]
                      << std::endl;
            exit(1);
        }
    }
    if (argc >= 7) {
        try {
            abs_mass_tol = boost::lexical_cast<double>(argv[6]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid abs_mass_tol (Expecting numerical): " << argv[6]
                      << std::endl;
            exit(1);
        }
    }
    if (argc >= 8) {
        output_filename = argv[7];
        to_stdout = false;
    }

    double cumulative_intensity_thresh = 100;
    double min_relative_intensity = 0.0;
    int min_peak_cutoffs = 1, clean_target_spectra = 0;
    int max_peak_cutoffs = 1000000;

    int quantise_spectra_dec_pl = -1;
    if (argc >= 9) {
        try {
            cumulative_intensity_thresh = boost::lexical_cast<double>(argv[8]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid cumulative_intensity_thresh (Expecting numerical): "
                      << argv[8] << std::endl;
            exit(1);
        }
    }
    if (argc >= 10) {
        try {
            min_relative_intensity = boost::lexical_cast<double>(argv[9]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid min_relative_intensity (Expecting numerical): "
                      << argv[9] << std::endl;
            exit(1);
        }
    }

    if (argc >= 11) {
        try {
            min_peak_cutoffs = boost::lexical_cast<int>(argv[10]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid min_peak_cutoffs (Expecting int) " << argv[10]
                      << std::endl;
            exit(1);
        }
    }
    if (argc >= 12) {
        try {
            max_peak_cutoffs = boost::lexical_cast<int>(argv[11]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid max_peak_cutoffs (Expecting int) " << argv[11]
                      << std::endl;
            exit(1);
        }
    }
    if (argc >= 13) {
        try {
            clean_target_spectra = boost::lexical_cast<bool>(argv[12]);
        } catch (boost::bad_lexical_cast &e) {
            std::cout << "Invalid apply_cutoffs (Expecting 0 or 1): " << argv[12]
                      << std::endl;
            exit(1);
        }
    }
    if (argc >= 14) {
        try {
            quantise_spectra_dec_pl = boost::lexical_cast<int>(argv[13]);
        } catch (boost::bad_lexical_cast e) {
            std::cout << "Invalid quantise_spectra_dec_pl (Expecting integer): "
                      << argv[13] << std::endl;
            exit(1);
        }
    }


    std::vector<MolData> data;
    parseInputFile(data, input_filename);

    // Prepare the score calculators
    Precision pcmp(ppm_mass_tol, abs_mass_tol);
    Recall rcmp(ppm_mass_tol, abs_mass_tol);
    WeightedRecall wrcmp(ppm_mass_tol, abs_mass_tol);
    WeightedPrecision wpcmp(ppm_mass_tol, abs_mass_tol);
    Jaccard jcmp(ppm_mass_tol, abs_mass_tol);
    WeightedJaccard wjcmp(ppm_mass_tol, abs_mass_tol);
    Dice dcmp(ppm_mass_tol, abs_mass_tol);
    WeightedDice wdcmp(ppm_mass_tol, abs_mass_tol);
    DotProduct dot(ppm_mass_tol, abs_mass_tol);
    OrigSteinDotProduct odot(ppm_mass_tol, abs_mass_tol);
    ScatterOutput shout(ppm_mass_tol, abs_mass_tol);

    // Set up the output (to file or stdout)
    std::streambuf *buf;
    std::ofstream of;
    if (!to_stdout) {
        of.open(output_filename.c_str());
        buf = of.rdbuf();
    } else
        buf = std::cout.rdbuf();
    std::ostream out(buf);

    // Out put configs
    out << "Config: " << std::endl;
    out << "cumulative_intensity_thresh " << cumulative_intensity_thresh << std::endl;
    out << "min_relative_intensity " << min_relative_intensity << std::endl;
    out << "peak_cutoffs(min, max) " << min_peak_cutoffs << " " << max_peak_cutoffs << std::endl;
    out << "clean_target_spectra " << clean_target_spectra << std::endl;
    out << "quantise_spectra_dec_pl " << quantise_spectra_dec_pl << std::endl << std::endl;


    // Compute the scores (per energy level and average scores across energy
    // levels)
    std::vector<double> rscores(data.size(), 0.0), pscores(data.size(), 0.0),
            jscores(data.size(), 0.0), wjscores(data.size(), 0.0),
            dscores(data.size(), 0.0), wdscores(data.size());

    std::vector<double> wrscores(data.size(), 0.0), wpscores(data.size(), 0.0),
            adscores(data.size(), 0.0);

    std::vector<double> dpscores(data.size(), 0.0), odpscores(data.size(), 0.0);

    std::vector<double> energy_rscores(data.size()), energy_pscores(data.size()),
            energy_jscores(data.size()), energy_dscores(data.size()),
            energy_wjscores(data.size()), energy_wdscores(data.size());

    std::vector<double> energy_wrscores(data.size()),
            energy_wpscores(data.size()), energy_adscores(data.size());
    std::vector<double> energy_dpscores(data.size()),
            energy_odpscores(data.size());

    double norm = 1.0 / (double) num_spectra;
    for (unsigned int i = 0; i < num_spectra; i++) {

        shout.out << "Energy" << i << std::endl;

        unsigned int idx = 0;
        std::vector<MolData>::iterator mit;
        for (mit = data.begin(); mit != data.end(); ++mit) {

            if (group_to_compute != -1 && mit->getGroup() != group_to_compute)
                continue;

            double mw = 0.0;
            try {
                mw = mit->getMolecularWeight();
            } catch (RDKit::MolSanitizeException e) {
                std::cout << "Could not sanitize mol " << mit->getId() << std::endl;
                continue;
            }

            // Read in the spectra (both measured and predicted)
            if (measured_is_msp)
                mit->readInSpectraFromMSP(*measured_msp);
            else {
                std::string spec_file = measured_spec_dir + "/" + mit->getId() + ".txt";
                mit->readInSpectraFromFile(spec_file);
            }
            if (predicted_is_msp) {
                try {
                    mit->readInSpectraFromMSP(*predicted_msp, true);
                } catch (MspIdException e) {
                    mit->freeSpectra();
                    continue;
                } // If it's not in the MSP, we failed to predict/enumerate it
            } else {
                std::string pred_spec_file =
                        predicted_spec_dir + "/" + mit->getId();
                if (boost::filesystem::exists( pred_spec_file + ".log" ) )
                    mit->readInSpectraFromFile(pred_spec_file + ".log" , true);
                else if (boost::filesystem::exists( pred_spec_file + ".txt" ) )
                    mit->readInSpectraFromFile(pred_spec_file + ".txt" , true);
                else
                    std::cerr << "Can not open " << pred_spec_file + ".txt" << " or " << pred_spec_file + ".log" << std::endl;

            }
            if (quantise_spectra_dec_pl >= 0) {
                mit->quantisePredictedSpectra(quantise_spectra_dec_pl);
                mit->quantiseMeasuredSpectra(quantise_spectra_dec_pl);
            }
            //if (apply_cutoffs)
            // mit->postprocessPredictedSpectra(cumulative_intensity_thresh, 5, 30, min_relative_intensity);
            //else
            // skip quantisePredictedSpectra since it is already doen
            mit->postprocessPredictedSpectra(cumulative_intensity_thresh, min_peak_cutoffs, max_peak_cutoffs,
                                             min_relative_intensity, -1);

            //num_spectra = mit->GetNumSpectra();
            if (clean_target_spectra)
                mit->cleanSpectra(0.1, 10.0);

            out << "\t" << mit->getSpectrum(i)->size() << "\t"
                << mit->getPredictedSpectrum(i)->size() << "\t";
            energy_rscores[idx] =
                    rcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            rscores[idx] += norm * energy_rscores[idx];
            energy_pscores[idx] =
                    pcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            pscores[idx] += norm * energy_pscores[idx];
            energy_wrscores[idx] =
                    wrcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            wrscores[idx] += norm * energy_wrscores[idx];
            energy_wpscores[idx] =
                    wpcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            wpscores[idx] += norm * energy_wpscores[idx];
            // jaccard Score
            energy_jscores[idx] =
                    jcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            jscores[idx] += norm * energy_jscores[idx];
            // weighted jaccard score
            energy_wjscores[idx] =
                    wjcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            wjscores[idx] += norm * energy_wjscores[idx];
            // dice coef
            energy_dscores[idx] =
                    dcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            dscores[idx] += norm * energy_dscores[idx];
            // weighted dice coef
            energy_wdscores[idx] =
                    wdcmp.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            wdscores[idx] += norm * energy_wdscores[idx];
            // dot product
            energy_dpscores[idx] =
                    dot.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            dpscores[idx] += norm * energy_dpscores[idx];
            energy_odpscores[idx] =
                    odot.computeScore(mit->getSpectrum(i), mit->getPredictedSpectrum(i));
            odpscores[idx] += norm * energy_odpscores[idx];
            out << mit->getId() << "\t" << energy_rscores[idx] << "\t"
                << energy_pscores[idx] << "\t" << energy_wrscores[idx] << "\t"
                << energy_wpscores[idx] << "\t" << energy_jscores[idx] << "\t" << energy_wjscores[idx] << "\t"
                << energy_dscores[idx] << "\t" << energy_wdscores[idx] << "\t"
                << energy_dpscores[idx] << "\t" << energy_odpscores[idx] << std::endl;
            idx++;
        }

        out << "Num mols in computation:" << idx << std::endl;
        energy_wpscores.resize(idx);
        energy_jscores.resize(idx);
        energy_wjscores.resize(idx);
        energy_dscores.resize(idx);
        energy_wdscores.resize(idx);
        energy_wrscores.resize(idx);
        energy_rscores.resize(idx);
        energy_pscores.resize(idx);
        energy_dpscores.resize(idx);
        energy_odpscores.resize(idx);
        wpscores.resize(idx);
        jscores.resize(idx);
        wjscores.resize(idx);
        dscores.resize(idx);
        wdscores.resize(idx);
        wrscores.resize(idx);
        rscores.resize(idx);
        pscores.resize(idx);
        adscores.resize(idx);
        dpscores.resize(idx);
        odpscores.resize(idx);

        // Report the mean and std of each energy
        out << "Energy " << i << std::endl;
        out << std::endl;
        out << "Recall (mean, std): ";
        reportMeanStd(out, energy_rscores);
        out << std::endl << "Precision (mean, std): ";
        reportMeanStd(out, energy_pscores);
        out << std::endl << "Weighted Recall (mean, std): ";
        reportMeanStd(out, energy_wrscores);
        out << std::endl << "Weighted Precision (mean, std): ";
        reportMeanStd(out, energy_wpscores);
        out << std::endl << "Jaccard (mean, std): ";
        reportMeanStd(out, energy_jscores);
        out << std::endl << "Weighted Jaccard (mean, std): ";
        reportMeanStd(out, energy_wjscores);
        out << std::endl << "Dice (mean, std): ";
        reportMeanStd(out, energy_dscores);
        out << std::endl << "Weighted Dice (mean, std): ";
        reportMeanStd(out, energy_wdscores);
        out << std::endl << "Dot Product (mean, std): ";
        reportMeanStd(out, energy_dpscores);
        out << std::endl << "Original Stein Dot Product (mean, std): ";
        reportMeanStd(out, energy_odpscores);
        out << std::endl;
    }


    // Write all the resulting scores to the output file (or stdout)
    out << "Totals:" << std::endl;

    // Report the mean and std of each
    out << std::endl;
    out << "Recall (mean, std): ";
    reportMeanStd(out, rscores);
    out << std::endl << "Precision (mean, std): ";
    reportMeanStd(out, pscores);
    out << std::endl << "Weighted Recall (mean, std): ";
    reportMeanStd(out, wrscores);
    out << std::endl << "Weighted Precision (mean, std): ";
    reportMeanStd(out, wpscores);
    out << std::endl << "Jaccard (mean, std): ";
    reportMeanStd(out, jscores);
    out << std::endl << "Weighted Jaccard (mean, std): ";
    reportMeanStd(out, wjscores);
    out << std::endl << "Dice (mean, std): ";
    reportMeanStd(out, dscores);
    out << std::endl << "Weighted Dice (mean, std): ";
    reportMeanStd(out, wdscores);
    //out << std::endl << "Altered Dot Product (mean, std): ";
    //reportMeanStd(out, adscores);
    out << std::endl << "Dot Product (mean, std): ";
    reportMeanStd(out, dpscores);
    out << std::endl << "Original Stein Dot Product (mean, std): ";
    reportMeanStd(out, odpscores);
    out << std::endl;

    if (measured_is_msp)
        delete measured_msp;
    if (predicted_is_msp)
        delete predicted_msp;

    return (0);
}

double addSq(double x, double y) { return x + y * y; }

void reportMeanStd(std::ostream &out, std::vector<double> &scores) {

    /*double mean =
            std::accumulate(scores.begin(), scores.end(), 0.0) / scores.size();
    double sq_mean =
            std::accumulate(scores.begin(), scores.end(), 0.0, addSq) / scores.size();
    double stddev = std::sqrt(sq_mean - mean * mean) / scores.size();*/

    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::variance> > acc;
    for (auto &score : scores)
        acc(score);

    out << mean(acc) << " " << sqrt(variance(acc)) << std::endl;
}

void parseInputFile(std::vector<MolData> &data, std::string &input_filename) {

    std::string line, smiles_or_inchi, id;
    std::ifstream ifs(input_filename.c_str(), std::ifstream::in);
    int group, num_mols = 0;
    config_t cfg;
    initDefaultConfig(cfg);

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
        data.push_back(MolData(id, smiles_or_inchi, group, &cfg));
    }
}
