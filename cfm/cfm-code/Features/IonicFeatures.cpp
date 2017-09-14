
#include "IonicFeatures.h"
void IonicFeatures::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
	
	int nl_pos = 0, nl_neg = 0, ion_pos = 0, ion_neg = 0;

	RDKit::ROMol::AtomIterator ai;
	for( ai = nl->mol.get()->beginAtoms(); ai != nl->mol.get()->endAtoms(); ++ai ){
		int ionic_frag_q; (*ai)->getProp("IonicFragmentCharge", ionic_frag_q);
		if( ionic_frag_q < 0 ) nl_neg = 1;
		if( ionic_frag_q > 0 ) nl_pos = 1;
	} 
	for( ai = ion->mol.get()->beginAtoms(); ai != ion->mol.get()->endAtoms(); ++ai ){
		int ionic_frag_q; (*ai)->getProp("IonicFragmentCharge", ionic_frag_q);
		if( ionic_frag_q < 0 ) ion_neg = 1;
		if( ionic_frag_q > 0 ) ion_pos = 1;
	} 
	
	fv.addFeature( nl_pos );	 	//NL has positive ionic fragment
	fv.addFeature( ion_pos );	//Ion has positive ionic fragment
	fv.addFeature( nl_neg );		//NL has negative ionic fragment
	fv.addFeature( ion_neg );	//Ion has negative ionic fragment
	fv.addFeature( !( nl_pos || nl_pos || ion_neg || ion_pos) );	//No ionic fragments anywhere
}
