void IonRootMMFFAtomType::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const
{
	int offset = fv.getTotalLength() - 1;
	int ring_break;
	nl->mol.get()->getProp( "IsRingBreak", ring_break );

	//Set features for atom types
	int atomtype, otheratomtype = 0;
	ion->root->getProp<int>("MMFFAtomType", atomtype);
	fv.addFeatureAtIdx(1.0, offset + atomtype );
	if( ring_break ){ 
		ion->other_root->getProp<int>("MMFFAtomType", otheratomtype);
		fv.addFeatureAtIdx(1.0, offset + otheratomtype );
	}
	//100 Features in total - last features indicates out-of-range
	if( atomtype < 1 || atomtype > 99 ) 
		fv.addFeatureAtIdx(1.0, offset + 100);
	else if(ring_break && ( otheratomtype < 1 || otheratomtype > 99 ) )
		fv.addFeatureAtIdx(1.0, offset + 100);
	else
		fv.addFeatureAtIdx(0.0, offset + 100);

}