#include "../Feature.h"

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
