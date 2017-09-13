void NLExtraFunctionalGroupFeatures::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	addFunctionalGroupFeatures( fv, nl, 1, ring_break, true );
}