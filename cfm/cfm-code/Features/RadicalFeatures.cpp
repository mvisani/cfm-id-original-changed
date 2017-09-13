void RadicalFeatures::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
    
        int ion_radical = moleculeHasSingleRadical( ion->mol.get() );
        int nl_radical = moleculeHasSingleRadical( nl->mol.get()  );
        fv.addFeature( ion_radical );						//Ion is radical
        fv.addFeature( nl_radical );						//NL is radical
        fv.addFeature( !ion_radical && !nl_radical );		//Neither NL or Ion are radical
    }
    