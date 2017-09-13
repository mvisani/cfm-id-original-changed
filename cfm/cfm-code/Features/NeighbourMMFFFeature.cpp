void NeighbourMMFFFeature::addNeighbourAtomTypes( FeatureVector &fv, const RootedROMolPtr *mol, const RDKit::Atom *root, int offset ) const
{
	//Iterate over the neighbours of the root atom
	RDKit::ROMol::ADJ_ITER_PAIR itp = mol->mol->getAtomNeighbors( root );
	int num_added = 0;
	for( ; itp.first != itp.second; ++itp.first ){

			RDKit::Atom *nbr_atom = mol->mol->getAtomWithIdx(*itp.first);
			int atomtype;
			nbr_atom->getProp<int>("MMFFAtomType", atomtype);
			fv.addFeatureAtIdx(1.0, offset + atomtype );
			if( atomtype < 1 || atomtype > 99 ) 
				fv.addFeatureAtIdx(1.0, offset + 100);
			num_added++;
	}
	//Additional feature indicating no neighbours
	if( num_added == 0 ) fv.addFeatureAtIdx(1.0, offset + 101);
}

