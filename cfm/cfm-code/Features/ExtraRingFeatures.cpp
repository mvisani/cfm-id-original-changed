#include "ExtraRingFeatures.h"

void ExtraRingFeatures::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const{
	
	//Not a ring break
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	fv.addFeature( !ring_break );	
	
	//Ion root is in ring
	RDKit::MolOps::findSSSR( *ion->mol );	
	RDKit::RingInfo *rinfo = ion->mol->getRingInfo();
	fv.addFeature( rinfo->minBondRingSize(ion->root->getIdx()) > 0 );

	//NL root is in ring
	RDKit::MolOps::findSSSR( *nl->mol );	
	rinfo = nl->mol->getRingInfo();
	fv.addFeature( rinfo->minBondRingSize(nl->root->getIdx()) > 0 );
}

