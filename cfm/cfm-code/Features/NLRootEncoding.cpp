void NLRootEncoding::compute( FeatureVector &fv, const RootedROMolPtr *nl ) const
{
    RDKit::ROMol &nl_ref = *(nl->mol.get());

    unsigned int minPath=1;
    unsigned int maxPath=7;
    unsigned int fpSize = 2048;
    ExplicitBitVect *fingerPrint = RDKit::RDKFingerprintMol(nl_ref, minPath, maxPath, fpSize);
    for(unsigned int i = 0; i < fingerPrint->getNumBits(); ++i)
    {
        fv.addFeature((*fingerPrint)[i]);
    }
}
