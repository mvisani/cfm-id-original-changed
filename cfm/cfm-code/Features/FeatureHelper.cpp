void FeatureHelper::initialiseRoots( RDKit::RWMol *rwmol ){
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai )
	{
		(*ai)->setProp("Root", 0 );
		(*ai)->setProp("OtherRoot", 0 );
	}
}

void FeatureHelper::labelGasteigers( RDKit::RWMol *rwmol ){
	// For each atom, will store the result in prop "OrigGasteigerCharge"
	RDKit::computeGasteigerCharges(rwmol); 
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai )
	{
		double gc;
		(*ai)->getProp<double>("_GasteigerCharge", gc);
		(*ai)->setProp("OrigGasteigerCharge", gc );
	}
}


void FeatureHelper::labelFunctionalGroups( RDKit::RWMol *rwmol, bool extra )
{

	const RDKit::MOL_SPTR_VECT &fgrps = fparams->getFuncGroups();
	const RDKit::MOL_SPTR_VECT &xfgrps = xfparams->getFuncGroups();

	std::vector<std::vector<unsigned int> > atom_fgidxs(rwmol->getNumAtoms());
	
	std::string prop_name;
	int num_grps, idx = 0;
    RDKit::MOL_SPTR_VECT::const_iterator fgrpi, fgrpe;
	if( extra ){ 
		fgrpi = xfgrps.begin(); fgrpe = xfgrps.end(); 
		num_grps = NUM_EXTRA_FGRPS; 
		prop_name = "ExtraFunctionalGroups";
	}
	else{ fgrpi = fgrps.begin(); fgrpe = fgrps.end(); 
		  num_grps = NUM_FGRPS; 
		prop_name = "FunctionalGroups";
	}

    for (; fgrpi != fgrpe; ++fgrpi, idx++) {
      std::string fname;
      (*fgrpi)->getProp("_Name", fname);
	  std::vector<RDKit::MatchVectType> fgpMatches;  //The format for each match is (queryAtomIdx, molAtomIdx)
      RDKit::SubstructMatch(*rwmol, *(fgrpi->get()), fgpMatches);
	  
	  std::vector<RDKit::MatchVectType>::const_iterator mat_it = fgpMatches.begin();
	  for (; mat_it != fgpMatches.end(); ++mat_it ) {
		  RDKit::MatchVectType::const_iterator it = (*mat_it).begin();
		  for( ; it != (*mat_it).end(); ++it )
			atom_fgidxs[it->second].push_back(idx);
	  }
	}
	// For each atom, store the list of functional group indexes in property "FunctionalGroups"
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai )
	{
		//Add an additional function group to indicate 'No Functional Groups'
		if( atom_fgidxs[(*ai)->getIdx()].size() == 0 )  
			atom_fgidxs[(*ai)->getIdx()].push_back( num_grps );
		(*ai)->setProp(prop_name, atom_fgidxs[(*ai)->getIdx()]);
	}

}

void FeatureHelper::labelMMFFAtomTypes( RDKit::RWMol *rwmol )
{
	
	//H-H causes exception...so assign to H-C atom type 5 (which is not used anyway)
	if( rwmol->getAtomWithIdx(0)->getSymbol() == "H" ){
		rwmol->getAtomWithIdx(0)->setProp("MMFFAtomType", (int)5 );
		return;
	}
	// For each atom, will store the result in prop "MMFFAtomType"
	RDKit::MMFF::MMFFMolProperties molprop( *rwmol );

	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){
		uint8_t atomtype = molprop.getMMFFAtomType( (*ai)->getIdx() );
		(*ai)->setProp("MMFFAtomType", (int)atomtype );
	}

}

void FeatureHelper::labelAromatics( RDKit::RWMol *rwmol ){	
	
	//Set bond aromaticity information
	RDKit::ROMol::BondIterator bi;
	for( bi = rwmol->beginBonds(); bi != rwmol->endBonds(); ++bi )
	{
		int aromatic = (*bi)->getIsAromatic();
		(*bi)->setProp("InAromaticRing", aromatic);
		(*bi)->setProp("InDblAromaticRing", 0);
	}

	//Check for any double-aromatic systems
	RDKit::MolOps::findSSSR( *rwmol );	
	RDKit::RingInfo *rinfo = rwmol->getRingInfo();
	std::vector<int> double_aromatic_idxs;
	for( unsigned int i = 0; i < rwmol->getNumBonds(); i++ )
	{
		if( rinfo->numBondRings(i) <= 1 ) continue;
		RDKit::Bond *bond = rwmol->getBondWithIdx(i);
		if( bond->getIsAromatic() ) 
			double_aromatic_idxs.push_back(i);
	}
	
	//If any are found, label all the bonds within them
	if(double_aromatic_idxs.size() == 0 ) return;
		
	//Consider each ring...
	RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
	RDKit::RingInfo::VECT_INT_VECT::iterator bit = brings.begin();
	for( ; bit != brings.end(); ++bit ){

		//Check for a double aromatic bond within the ring
		bool hasDblArom = false;
		RDKit::RingInfo::INT_VECT::iterator it;
		for( it = bit->begin(); it != bit->end(); ++it ){
			std::vector<int>::iterator ii = double_aromatic_idxs.begin();
			for( ; ii != double_aromatic_idxs.end(); ++ii )
				if( *ii == *it ) hasDblArom = true;
			if(hasDblArom) break;
		}

		//If one exists, label all bonds in the ring
		if( !hasDblArom ) continue;
		for( it = bit->begin(); it != bit->end(); ++it ){
			RDKit::Bond *bond = rwmol->getBondWithIdx(*it);
			bond->setProp("InDblAromaticRing", 1);
		}
		
	}
}

void FeatureHelper::labelOriginalMasses( RDKit::RWMol *rwmol )
{
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai )
	{
		double mass = 0.0;
		std::string symbol = (*ai)->getSymbol();
		mass += pt->getMostCommonIsotopeMass(symbol);
		mass += (*ai)->getTotalNumHs()*pt->getMostCommonIsotopeMass("H");
		(*ai)->setProp("OriginalMass", mass);
	}
}


void FeatureHelper::labelAtomsWithLonePairs( RDKit::RWMol *rwmol )
{
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
	RDKit::ROMol::AtomIterator ai;
	RDKit::MolOps::findSSSR( *rwmol );	
	RDKit::RingInfo *rinfo = rwmol->getRingInfo();
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){
		std::string symbol = (*ai)->getSymbol();
		int nouter = pt->getNouterElecs( symbol.c_str() );
		int def_val = pt->getDefaultValence(symbol.c_str());
		(*ai)->setProp("HasLP", (int)(nouter > def_val && def_val != 1 && def_val != -1) );	//Allow O,N,S,P..but not C, Halogens, Metals,.
	}
}

void FeatureHelper::labelOriginalBondTypes( RDKit::RWMol *rwmol )
{	
	RDKit::ROMol::BondIterator bi;	
	for( bi = rwmol->beginBonds(); bi != rwmol->endBonds(); ++bi ){
		if( (*bi)->getIsAromatic() ) (*bi)->setProp("OrigBondType", 4);
		else if( (*bi)->getIsConjugated() ) (*bi)->setProp("OrigBondType", 5);
		else (*bi)->setProp("OrigBondType", (int)((*bi)->getBondTypeAsDouble()) );
	}
}