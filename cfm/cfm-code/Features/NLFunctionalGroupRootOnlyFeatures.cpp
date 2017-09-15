void NLFunctionalGroupRootOnlyFeatures::compute(
    FeatureVector &fv, const RootedROMolPtr *ion,
    const RootedROMolPtr *nl) const {
  int ring_break;
  nl->mol.get()->getProp("IsRingBreak", ring_break);
  addFunctionalGroupFeatures(fv, nl, 0, ring_break);
}