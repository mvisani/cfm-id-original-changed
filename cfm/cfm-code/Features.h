/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Features.h
#
# Description: 	Code for computing features for fragmentations.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FEATURE_H__
#define __FEATURE_H__

#include "Util.h"
#include "FunctionalGroups.h"

#include <GraphMol/FragCatalog/FragCatParams.h>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/lexical_cast.hpp>

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>


struct input_file_t;

//Exception to throw when the input feature configuration file is invalid 
class InvalidConfigException: public std::exception{

	virtual const char* what() const throw(){
		return "Invalid Feature Configuration File";
	}
};

class FeatureCalculationException: public std::exception{
private:
    std::string message_;
public:
	FeatureCalculationException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << "Error computing feature vector: " << message_ << std::endl;
		return message_.c_str();
	}
	~FeatureCalculationException() throw() {};
};

class FeatureHelperException: public std::exception{
private:
    std::string message_;
public:
	FeatureHelperException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << "Error in FeatureHelper: " << message_ << std::endl;
		return message_.c_str();
	}
	~FeatureHelperException() throw() {};
};

//Structure to hold a sparse computed feature vector
typedef unsigned int feature_t;
class FeatureVector{
public:
	FeatureVector(){fv_idx = 0;};
	void addFeature( double value );
	void addFeatureAtIdx( double value, unsigned int idx );
	unsigned int const getTotalLength() const {return fv_idx;};
	feature_t const getFeature( int idx ) const { return fv[idx];};
	std::vector<feature_t>::const_iterator getFeatureBegin() const { return fv.begin(); };
	std::vector<feature_t>::const_iterator getFeatureEnd() const { return fv.end(); };
	unsigned int const getNumSetFeatures() const { return fv.size();};

private:
	std::vector<feature_t> fv;
	unsigned int fv_idx;
};

//Base class to compute a feature - all features should inherit from this
class Feature{

public:
	virtual void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const = 0;
	unsigned int getSize() const { return size; };
	std::string getName() const { return name; };
	virtual ~Feature(){};

protected:
	unsigned int size;	
	std::string name;
	static const std::vector<std::string> &OKsymbols();
	static const std::vector<std::string> &OKSymbolsLess();
	void replaceUncommonWithX( std::string &symbol ) const;
	int getSymbolsLessIndex( std::string& symbol) const;

};

//Class to compute a feature vector
class FeatureCalculator{

public:
	//Constructor: Initialise the calculator using a config file listing features
	FeatureCalculator( std::string &config_filename );

	//Constructor: Initialise the calculator using a list of feature names
	FeatureCalculator( std::vector<std::string> &feature_list );
	
	//Compute the expected number of total features
	unsigned int getNumFeatures();
	
	//Retrieve the list of feature names being used
	std::vector<std::string> getFeatureNames();

	//Retrieve a list of valid feature names (for testing)
	static const std::vector<std::string> getValidFeatureNames();

	//Compute the feature vector for the input ion and nl (with labeled Root atoms)
	// - NB: responsibility of caller to delete.
	FeatureVector *computeFV( const RootedROMolPtr *ion, const RootedROMolPtr *nl );

	bool includesFeature( const std::string &fname );

private:
	//List of feature classes ready to be used
	static const boost::ptr_vector<Feature> &featureCogs();

	//Indexes of feature classes that are selected for use
	std::vector<int> used_feature_idxs;

	//Helper function - Configure feature for use
	void configureFeature( std::string &name );
};

//*************************
//FEATURE IMPLEMENTATIONS:
//*************************

class BreakAtomPair : public Feature {
public:
	BreakAtomPair(){ size = 72; name = "BreakAtomPair"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class RootPathFeature : public Feature {
protected:
	typedef std::vector<std::string> path_t;
	// function to compute path with given length from a root
	void computeRootPaths(std::vector<path_t> &paths, 
							const RootedROMolPtr *mol, 
							int len, 
							bool ring_break,
							bool with_bond) const;

	// function to add features with a length of two						
	void addRootPairFeatures(FeatureVector &fv, 
							std::vector<path_t> &paths, 
							int ring_break) const;

	// function to add features with a length of three
	void addRootTripleFeatures(FeatureVector &fv, 
								std::vector<path_t> &paths, 
								int ring_break,
								bool with_bond) const;
private:
	// function to add path from given atom
	void addPathsFromAtom( std::vector<path_t> &paths, const RDKit::Atom *atom, 
						   const romol_ptr_t mol, 
						   const RDKit::Atom *prev_atom, 
						   path_t &path_so_far, 
						   int len, 
						   bool with_bond) const;
};



class IonRootTriples : public RootPathFeature {
public:
	IonRootTriples(){ size = 865; name = "IonRootTriples"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootTriples : public RootPathFeature {
public:
	NLRootTriples(){ size = 865; name = "NLRootTriples"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class IonRootTriplesIncludeBond : public RootPathFeature {
public:
	IonRootTriplesIncludeBond
	(){ size = 10585; name = "IonRootTriplesIncludeBond"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootTriplesIncludeBond : public RootPathFeature {
public:
	NLRootTriplesIncludeBond(){ size = 10585; name = "NLRootTriplesIncludeBond"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};



// Features use fingerprint encode ion fragmentation
class IonRootEncoding : public RootPathFeature {
	IonRootEncoding(){ size = 2048; name = "IonRootEncoding"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion ) const;
};

// Features use fingerprint encode Neutral Loss
class NLRootEncoding : public RootPathFeature {
	NLRootEncoding(){ size = 2048; name = "NLRootEncoding"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *nl ) const;
};



class BrokenOrigBondType : public Feature {
public:
	BrokenOrigBondType(){ size = 7; name = "BrokenOrigBondType"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NeighbourOrigBondTypes : public Feature {
public:
	NeighbourOrigBondTypes(){ size = 12; name = "NeighbourOrigBondTypes"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class IonicFeatures: public Feature {
public:
	IonicFeatures(){ size = 5; name = "IonicFeatures"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};


class ExtraRingFeatures : public Feature {
public:
	ExtraRingFeatures(){ size = 3; name = "ExtraRingFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class IonRootMMFFAtomType : public Feature {
public:
	IonRootMMFFAtomType(){ size = 100; name = "IonRootMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootMMFFAtomType : public Feature {
public:
	NLRootMMFFAtomType(){ size = 100; name = "NLRootMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};



class IonRootAtom : public RootAtomFeature {
public:
	IonRootAtom(){ size = 13; name = "IonRootAtom";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootAtom : public RootAtomFeature {
public:
	NLRootAtom(){ size = 13; name = "NLRootAtom";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};


class NeighbourMMFFFeature : public Feature {
protected:
	void addNeighbourAtomTypes( FeatureVector &fv, const RootedROMolPtr *mol, const RDKit::Atom *root, int offset ) const;

};

class IonNeighbourMMFFAtomType : public NeighbourMMFFFeature {
public:
	IonNeighbourMMFFAtomType(){ size = 101; name = "IonNeighbourMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLNeighbourMMFFAtomType : public NeighbourMMFFFeature {
public:
	NLNeighbourMMFFAtomType(){ size = 101; name = "NLNeighbourMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};



class NLFunctionalGroupFeaturesD2 : public FunctionalGroupFeature {
public:
	NLFunctionalGroupFeaturesD2(){ size = (NUM_FGRPS+1)*3; name = "NLFunctionalGroupFeaturesD2";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};


class NLFunctionalGroupFeatures : public FunctionalGroupFeature {
public:
	NLFunctionalGroupFeatures(){ size = (NUM_FGRPS+1)*2; name = "NLFunctionalGroupFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};


class NLFunctionalGroupRootOnlyFeatures : public FunctionalGroupFeature {
public:
	NLFunctionalGroupRootOnlyFeatures(){ size = NUM_FGRPS+1; name = "NLFunctionalGroupRootOnlyFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};



class NLExtraFunctionalGroupFeatures : public FunctionalGroupFeature {
public:
	NLExtraFunctionalGroupFeatures(){ size = (NUM_EXTRA_FGRPS+1)*2; name = "NLExtraFunctionalGroupFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};


class QuadraticFeatures : public Feature {
public:
	QuadraticFeatures(){ size = 0; name = "QuadraticFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};




#endif // __FEATURE_H__