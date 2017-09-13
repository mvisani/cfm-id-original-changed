void IonRootAtom::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
    
        int ring_break;
        nl->mol.get()->getProp( "IsRingBreak", ring_break );
        computeRootAtomFeature( fv, ion, ring_break );
    }
    