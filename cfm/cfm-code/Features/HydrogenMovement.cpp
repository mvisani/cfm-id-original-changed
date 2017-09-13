void HydrogenMovement::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl  ) const{
    
        double h_movement = 0.0;
    
        //Compute the mass difference in the ion
        RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();
        RDKit::ROMol::AtomIterator ai;
        for( ai = ion->mol.get()->beginAtoms(); ai != ion->mol.get()->endAtoms(); ++ai ){
            double orig_mass, mass = 0.0;
            std::string symbol = (*ai)->getSymbol();
            mass += pt->getMostCommonIsotopeMass(symbol);
            mass += (*ai)->getTotalNumHs()*pt->getMostCommonIsotopeMass("H");
            if( !(*ai)->hasProp("OriginalMass") ) std::cout << "No OriginalMass prop..." << std::endl;
            (*ai)->getProp<double>("OriginalMass", orig_mass);
            h_movement += (mass - orig_mass);
        }
    
        //Binary on/off indicating whether a particular transfer occurred
        for( double h = -4.0; h <= 4.0; h += 1.0 ){	
            if( fabs( h - h_movement ) < 0.5 ) fv.addFeature( 1.0 );
            else fv.addFeature(0.0);
        }
        //Catch-all for all other hydrogen transfers
        if( fabs( h_movement ) > 4.0 ) fv.addFeature( 1.0 );
        else fv.addFeature( 0.0 );
    }
    
    
    