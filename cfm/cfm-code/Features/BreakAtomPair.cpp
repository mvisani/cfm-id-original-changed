
#include "BreakAtomPair.h"

void BreakAtomPair::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
	
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	std::vector<symbol_pair_t> pairs;

	//Ion Symbol(s)
	std::string irootsymbol, iotherrootsymbol;
	irootsymbol = ion->root->getSymbol();
	replaceUncommonWithX(irootsymbol);
	if( ring_break ){
		iotherrootsymbol = ion->other_root->getSymbol();
		replaceUncommonWithX(iotherrootsymbol);
	}

	//Neutral Loss Symbol(s)
	std::string nlrootsymbol, nlotherrootsymbol;
	nlrootsymbol = nl->root->getSymbol();
	replaceUncommonWithX(nlrootsymbol);
	if( ring_break ){
		nlotherrootsymbol = nl->other_root->getSymbol();
		replaceUncommonWithX(nlotherrootsymbol);
	}

	//Pairs
	pairs.push_back( symbol_pair_t( irootsymbol, nlrootsymbol ));
	if( ring_break ) 
		pairs.push_back( symbol_pair_t( iotherrootsymbol, nlotherrootsymbol ));

	//Iterate through all combinations of atom pairs, appending
	//a feature for each; 1 if it matches, 0 otherwise.
	//Note: the order matters here, ion first then nl
	std::vector<std::string>::const_iterator it1, it2;
	const std::vector<std::string> *ok_symbols = &OKSymbolsLess();
	for( it1 = ok_symbols->begin(); it1 != ok_symbols->end(); ++it1 ){
		for( it2 = ok_symbols->begin(); it2 != ok_symbols->end(); ++it2 ){
			symbol_pair_t sp = symbol_pair_t(*it1, *it2);
			double nonringf = 0.0, ringf = 0.0;
			if( sp == *pairs.begin() ){ 
				nonringf = !ring_break;
				ringf = ring_break;
			}
			if( sp == *pairs.rbegin() ) ringf = ring_break;
			fv.addFeature(nonringf);
			fv.addFeature(ringf);
		} 
	}

}