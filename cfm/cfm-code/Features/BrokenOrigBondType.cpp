#include "BrokenOrigBondType.h"

void BrokenOrigBondType::compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const{
    
        int bondtype; nl->mol.get()->getProp( "BrokenOrigBondType", bondtype );
        fv.addFeature( bondtype == 1 );				//SINGLE
        fv.addFeature( bondtype == 2 );				//DOUBLE
        fv.addFeature( bondtype == 3 );				//TRIPLE
        fv.addFeature( bondtype == 4 );				//AROMATIC
        fv.addFeature( bondtype == 5 );				//CONJUGATED
        fv.addFeature( bondtype == 6 );				//IONIC
        fv.addFeature( bondtype == 7 );				//H ONLY
    
    }