/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Spectrum.cpp
#
# Description: 	Class for spectrum.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Spectrum.h"
#include "Util.h"
#include "Version.h"

#include <algorithm>
#include <cmath>

void Spectrum::outputToStream(std::ostream &out, bool do_annotate, int mz_precision, bool normalize_to_max) const {

	double max_intensity = normalize_to_max ? getMaxIntensity() : -1.0;

	for (auto &peak : peaks) {

		// compute display_intensity , and display if intensity is large enough
		double display_intensity_value = normalize_to_max ? peak.intensity / max_intensity * 100.0 : peak.intensity;
		int display_intensity          = std::floor(display_intensity_value * 100 + 0.5);
		if (display_intensity > 0) {
			out << std::fixed << std::setprecision(mz_precision) << peak.mass << " " << std::setprecision(2)
			    << display_intensity / 100.0;

			if (do_annotate) {
				std::stringstream ss_values;
				ss_values << std::setprecision(5) << "(";
				auto ita = peak.annotations.begin();
				for (; ita != peak.annotations.end(); ++ita) {
					out << " " << ita->first;
					if (ita != peak.annotations.begin()) ss_values << " ";
					ss_values << ita->second * 100.0;
				}
				ss_values << ")";
				if (!peak.annotations.empty()) out << " " << ss_values.str();
			}

			out << std::endl;
		}
	}
}

void Spectrum::getDisplayedFragmentIds(std::set<int> &ids, bool normalize_to_max) const {

	double max_intensity = normalize_to_max ? getMaxIntensity() : -1.0;

	for (auto &peak : peaks) {

		// compute display_intensity , and display if intensity is large enough
		double display_intensity_value = normalize_to_max ? peak.intensity / max_intensity * 100.0 : peak.intensity;
		int display_intensity          = std::floor(display_intensity_value * 100 + 0.5);
		if (display_intensity > 0) {
			for (auto &annotation : peak.annotations) { ids.insert(annotation.first); }
		}
	}
}
void Spectrum::outputToMspStream(std::ostream &out, std::string id, int ionization_mode, int energy,
                                 std::string &smiles_or_inchi, int mz_precision) const {

	if (ionization_mode == POSITIVE_EI_IONIZATION_MODE)
		out << "Name: +ve in-silico MS by ";
	else if (ionization_mode == POSITIVE_ESI_IONIZATION_MODE)
		out << "Name: +ve in-silico MS/MS by ";
	else
		out << "Name: -ve in-silico MS/MS by ";
	out << APP_STRING << " " << PROJECT_VER << " for " << id << std::endl;
	out << "ID: " << id << std::endl;
	out << "Smiles/Inchi:" << smiles_or_inchi << std::endl;
	out << "Comment: Energy" << energy << std::endl;
	out << "Num peaks: " << peaks.size() << std::endl;
	outputToStream(out, false, mz_precision, true);
	out << std::endl;
}

void Spectrum::outputToMgfStream(std::ostream &out, std::string id, int ionization_mode, int energy, double mw,
                                 std::string &smiles_or_inchi, int mz_precision) const {

	out << "BEGIN IONS" << std::endl;
	out << "PEPMASS=" << std::setprecision(10) << mw << std::endl;
	if (ionization_mode == POSITIVE_ESI_IONIZATION_MODE || ionization_mode == POSITIVE_EI_IONIZATION_MODE)
		out << "CHARGE=1+" << std::endl;
	else if (ionization_mode == NEGATIVE_ESI_IONIZATION_MODE)
		out << "CHARGE=1-" << std::endl;
	out << "TITLE=" << id << ";Energy" << energy << ";";
	if (ionization_mode == POSITIVE_EI_IONIZATION_MODE)
		out << "[M]+;In-silico MS by ";
	else if (ionization_mode == POSITIVE_ESI_IONIZATION_MODE)
		out << "[M+H]+In-silico MS/MS by ";
	else if (ionization_mode == NEGATIVE_ESI_IONIZATION_MODE)
		out << "[M-H]-;In-silico MS/MS by ";
	out << APP_STRING << " " << PROJECT_VER << ";" << smiles_or_inchi << ";" << std::endl;
	outputToStream(out, false, mz_precision, true);
	out << "END IONS" << std::endl;
}

void Spectrum::quantisePeaksByMass(int num_dec_places) {

	// Combine peaks that have the same mass when reduced to num_dec_places
	// Note: this is mostly used for extreme cases like the NIST data, where
	// masses are given only to integer precision.
	normalizeAndSort();

	long long prev_mass = 0;
	auto it             = peaks.begin();
	for (; it != peaks.end(); ++it) {
		auto tmp_mass = (long long)(lround(it->mass * std::pow(10.0, num_dec_places)));
		it->mass      = tmp_mass * std::pow(10.0, -num_dec_places);
		if (tmp_mass == prev_mass && it != peaks.begin()) {
			it->intensity += (it - 1)->intensity;
			it->annotations.insert(it->annotations.end(), (it - 1)->annotations.begin(), (it - 1)->annotations.end());
			it = peaks.erase(it - 1);
		}
		prev_mass = tmp_mass;
	}
	normalizeAndSort();
}

void Spectrum::postProcess(double perc_thresh, int min_peaks, int max_peaks, double min_relative_intensity_prec) {

	if (peaks.empty()) return;

	std::sort(peaks.begin(), peaks.end(), sort_peaks_by_intensity);
	auto max_intensity = peaks[0].intensity;
	double total       = 0.0;
	int count          = 0;
	// we need at least 1 peak
	min_peaks          = std::max(min_peaks, 1);
	for (auto &peak : peaks) {
		total += peak.intensity;
		// e.g. Take the top 80% of energy (assuming at least 5 peaks),
		// or the highest 30 peaks (whichever comes first)
		if ((total > perc_thresh && count > min_peaks) || count > max_peaks ||
		    (peak.intensity / max_intensity) * 100.0 < min_relative_intensity_prec) {
			break;
		}
		count++;
	}
	peaks.resize(count);
	std::sort(peaks.begin(), peaks.end(), sort_peaks_by_mass);
}

void Spectrum::normalizeAndSort() {

	if (!is_normalized) {
		// Compute the normalizer
		double sum = 0.0;
		auto itp   = peaks.begin();
		for (; itp != peaks.end(); ++itp) sum += itp->intensity;
		double norm = 1.0;
		if (sum > 0.0) norm = 100.0 / sum;

		// Adjust the values
		for (itp = peaks.begin(); itp != peaks.end(); ++itp) itp->intensity *= norm;
	}

	// Ensure the peaks are sorted by mass
	if (!is_sorted) std::sort(peaks.begin(), peaks.end(), sort_peaks_by_mass);

	is_sorted     = true;
	is_normalized = true;
}

void Spectrum::clean(double abs_mass_tol, double ppm_mass_tol) {

	// Ensure the initial spectrum is normalized and sorted by mass
	normalizeAndSort();

	// Filter peaks that are too close together
	// - Remove peaks within abs_mass_tol of a larger peak
	//	 and any peaks below an absolute intensity threshold
	double abs_intensity_thresh = 0.01;
	std::vector<bool> peak_flags(peaks.size(), true);

	// Forward Pass
	double prev_intensity = -1.0;
	double prev_mass      = -100.0;
	auto itp              = peaks.begin();
	for (int idx = 0; itp != peaks.end(); ++itp, idx++) {
		if (itp->intensity < abs_intensity_thresh) peak_flags[idx] = false;

		double mass_tol = getMassTol(abs_mass_tol, ppm_mass_tol, itp->mass);
		while (fabs(itp->mass - prev_mass) < mass_tol && itp->intensity < prev_intensity) {
			peak_flags[idx] = false;
			if (idx < peaks.size() - 1) {
				idx += 1;
				itp += 1;
			}
		}
		idx -= 1;
		itp -= 1;

		prev_mass      = itp->mass;
		prev_intensity = itp->intensity;
	}

	// Reverse Pass
	prev_intensity = -1.0;
	prev_mass      = -100.0;
	auto ritp      = peaks.rbegin();
	for (int idx = peaks.size() - 1; ritp != peaks.rend(); ++ritp, idx--) {
		double mass_tol = getMassTol(abs_mass_tol, ppm_mass_tol, ritp->mass);
		while (fabs(ritp->mass - prev_mass) < mass_tol && ritp->intensity < prev_intensity) {
			peak_flags[idx] = false;
			if (idx > 0) {
				idx -= 1;
				ritp += 1;
			}
		}
		idx += 1;
		ritp -= 1;

		prev_mass      = ritp->mass;
		prev_intensity = ritp->intensity;
	}

	// Alter the spectrum to include the selected peaks
	std::vector<Peak> peaks_copy(peaks);
	peaks.clear();
	itp = peaks_copy.begin();
	for (int idx = 0; itp != peaks_copy.end(); ++itp, idx++)
		if (peak_flags[idx]) peaks.push_back(*itp);

	// Re-normalize
	normalizeAndSort();
}

void Spectrum::sortAndNormalizeAnnotations() {

	auto it = peaks.begin();
	for (; it != peaks.end(); ++it) {
		// Sort
		std::sort(it->annotations.begin(), it->annotations.end(), sort_annotations_by_score);

		// Normalize
		auto itt     = it->annotations.begin();
		double total = 0.0;
		for (; itt != it->annotations.end(); ++itt) total += itt->second;
		double norm = 1.0;
		if (total > 0.0) norm = it->intensity / (100.0 * total);
		for (itt = it->annotations.begin(); itt != it->annotations.end(); ++itt) itt->second *= norm;
	}
}

int Spectrum::removePeaksWithNoFragment(std::vector<double> &frag_masses, double abs_tol, double ppm_tol) {

	int num_removed = 0;
	// Remove any peaks more than mass_tol away from any fragment
	auto peak       = peaks.begin();
	for (; peak != peaks.end();) {

		double mass_tol = getMassTol(abs_tol, ppm_tol, peak->mass);
		bool found      = false;
		for (auto frag_mass : frag_masses) {
			if (fabs(frag_mass - peak->mass) < mass_tol) {
				found = true;
				break;
			}
		}
		if (!found) {
			num_removed++;
			peak = peaks.erase(peak);
		} else
			++peak;
	}
	// Renormalise
	normalizeAndSort();

	return num_removed;
}

void Spectrum::convertToLogScale() {
	// find max intensity
	auto max_intensity = getMaxIntensity();

	for (auto &peak : peaks) max_intensity = std::max(peak.intensity, max_intensity);

	// convert peak height to relative to main peak (height 100.0)
	// make sure we don't have negative value
	for (auto &peak : peaks) peak.intensity = log(std::max(peak.intensity / max_intensity * 100.0, 1.0));
	normalizeAndSort();
}

void Spectrum::convertToLinearScale() {

	// find max intensity
	auto max_intensity = getMaxIntensity();

	auto ratio = 100.0 / log(101);
	for (auto &peak : peaks) peak.intensity = exp((peak.intensity / max_intensity) / ratio) - 1;
	normalizeAndSort();
}

double Spectrum::getMaxIntensity() const {
	auto max_intensity = 0.0;
	for (auto &peak : peaks) max_intensity = std::max(peak.intensity, max_intensity);
	return max_intensity;
}
