void IonNeighbourMMFFAtomType::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const
{
	int offset = fv.getTotalLength() - 1;
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );
	fv.addFeatureAtIdx(0.0, offset + 101);	//Make the feature vector the right length
	addNeighbourAtomTypes( fv, ion, ion->root, offset );
	if( ring_break ) addNeighbourAtomTypes( fv, ion, ion->other_root, offset );
	
}
