void NeighbourOrigBondTypes::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
    
        int ring_break;
        nl->mol.get()->getProp( "IsRingBreak", ring_break );
        addNeighbourOrigBondFeatures( fv, ion, ring_break );
        addNeighbourOrigBondFeatures( fv, nl, ring_break );
    }
    