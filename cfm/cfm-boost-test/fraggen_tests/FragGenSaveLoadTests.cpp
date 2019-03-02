/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for features.cpp
#
# Author: Felicity Allen, Fei Wang
# Created: November 2012, 2019
#########################################################################*/
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
namespace bdata = boost::unit_test::data;

#include "FragGenTestsUtils.h"

BOOST_AUTO_TEST_SUITE(FragGenSaveAndLoadTests)


    std::vector<std::string> smiles_or_inchis = {"CC1=CC(=O)CC(C)C1", "NCCCN"};

    BOOST_DATA_TEST_CASE(NoIsotopesSaveAndLoad, bdata::make(smiles_or_inchis), smiles_or_inchi){
        double tol = 1e-6;
        config_t cfg; initDefaultConfig( cfg );
        cfg.include_isotopes = false;
        MolData morig("Test ID", smiles_or_inchi.c_str(), &cfg );
        std::vector<std::string> fnames;
        fnames.push_back( "BreakAtomPair" );
        fnames.push_back( "HydrogenRemoval" );
        FeatureCalculator fc( fnames );
        morig.computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);

        //Save state
        std::string fvfilename = "tmp_fv_file.fg";
        if( boost::filesystem::exists( fvfilename ) )
            boost::filesystem::remove( fvfilename );
        morig.writeFVFragmentGraph(fvfilename);

        //Load state into another moldata object
        MolData mload("Test ID", smiles_or_inchi.c_str(), &cfg );
        mload.readInFVFragmentGraph(fvfilename);

        //Compare orig vs loaded
        BOOST_CHECK_EQUAL( morig.getNumFragments(),mload.getNumFragments());
        BOOST_CHECK_EQUAL( morig.getNumTransitions(),mload.getNumTransitions());
        BOOST_CHECK_EQUAL( morig.getFGHeight(), mload.getFGHeight());

        for( int i = 0; i < morig.getNumFragments(); i++ ){
            const Fragment *f1 = morig.getFragmentAtIdx(i);
            const Fragment *f2 = mload.getFragmentAtIdx(i);

            // Check id
            BOOST_CHECK_EQUAL( f1->getId(),f2->getId());
            // check mass
            BOOST_CHECK_CLOSE_FRACTION(f1->getMass(),  f2->getMass(), tol);
            // check flag
            BOOST_CHECK_EQUAL( f1->isIntermediate(),f2->isIntermediate());
        }
        for( int i = 0; i < morig.getNumTransitions(); i++ ){
            const TransitionPtr t1 = morig.getTransitionAtIdx(i);
            const TransitionPtr t2 = mload.getTransitionAtIdx(i);

            // Check from id
            BOOST_CHECK_EQUAL( t1->getFromId(), t2->getFromId());
            // check to id
            BOOST_CHECK_EQUAL( t1->getToId(), t2->getToId());
            // check flag

            const FeatureVector *fv1 = morig.getFeatureVectorForIdx(i);
            const FeatureVector *fv2 = mload.getFeatureVectorForIdx(i);
            // check fv size
            BOOST_CHECK_EQUAL( fv1->getTotalLength(), fv2->getTotalLength());
            BOOST_CHECK_EQUAL( fv1->getNumSetFeatures(), fv2->getNumSetFeatures());

            auto fit1 = fv1->getFeatureBegin();
            auto fit2 = fv2->getFeatureBegin();
            for(; fit1 != fv1->getFeatureEnd(); ++fit1, ++fit2 ){
                BOOST_CHECK_EQUAL( *fit1, *fit2 );
            }
        }

        //Check the tmaps
        const tmap_t *fromidmap1 = morig.getFromIdTMap();
        const tmap_t *fromidmap2 = mload.getFromIdTMap();
        const tmap_t *toidmap1 = morig.getToIdTMap();
        const tmap_t *toidmap2 = mload.getToIdTMap();

        BOOST_CHECK_EQUAL( fromidmap1->size(), fromidmap2->size());
        BOOST_CHECK_EQUAL( toidmap1->size(), toidmap2->size());

        for( int i = 0; i < fromidmap1->size(); i++ ){
            BOOST_CHECK_EQUAL_COLLECTIONS((*fromidmap1)[i].begin(), (*fromidmap1)[i].end(),
                                          (*fromidmap2)[i].begin(), (*fromidmap2)[i].end());

            BOOST_CHECK_EQUAL_COLLECTIONS((*toidmap1)[i].begin(), (*toidmap1)[i].end(),
                                          (*toidmap2)[i].begin(), (*toidmap2)[i].end());
        }
    }

BOOST_AUTO_TEST_SUITE_END()