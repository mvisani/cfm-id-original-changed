#include "../Features.h"

class BreakAtomPair : public Feature {
public:
	BreakAtomPair(){ size = 72; name = "BreakAtomPair"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};