void IonRootTriplesIncludeBond::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
	
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	std::vector<path_t> paths;
	computeRootPaths( paths, ion, 3, ring_break, true);
	addRootTripleFeatures( fv, paths, ring_break, true);
}