void NLRootPairs::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
    
        int ring_break;
        nl->mol.get()->getProp( "IsRingBreak", ring_break );
        std::vector<path_t> paths;
        computeRootPaths( paths, nl, 2, ring_break, false);
        addRootPairFeatures( fv, paths, ring_break);
    }
    