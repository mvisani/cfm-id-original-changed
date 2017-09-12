void IonRootEncoding::compute( FeatureVector &fv, const RootedROMolPtr *ion ) const
{
    RDKit::ROMol &ion_ref = *(ion->mol.get());

    unsigned int minPath=1;
    unsigned int maxPath=7;
    unsigned int fpSize = 2048;
    ExplicitBitVect *fingerPrint = RDKit::RDKFingerprintMol(ion_ref, minPath, maxPath, fpSize);
    for(unsigned int i = 0; i < fingerPrint->getNumBits(); ++i)
    {
        fv.addFeature((*fingerPrint)[i]);
    }
}