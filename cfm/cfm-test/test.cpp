/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# test.cpp
#
# Description: 	Test code.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "mpi.h"
#include "test.h"

#include <iostream> 
#include <exception>
#include <vector>

#include "tests/message_test.h"
#include "tests/features_test.h"
#include "tests/fraggen_test.h"
#include "tests/param_test.h"
#include "tests/nn_param_test.h"
#include "tests/inference_test.h"
#include "tests/libdai_inference_test.h"
#include "tests/comms_test.h"
#include "tests/em_test.h"
#include "tests/comparator_test.h"
#include "tests/spectrum_test.h"
#include "tests/msp_reader_test.h"

#include <boost/ptr_container/ptr_vector.hpp>

Test::Test(){
	description = "None";
}

void Test::run(){
	try{
		excpt_occurred = false;
		passed = false;
		runTest();
	}
	catch(std::exception & e){
		excpt_occurred = true;
		std::cout << e.what();
		throw e;
	}
}

void Test::runTest(){
	std::cout << "Running parent runTest" << std::endl;
	//Do nothing: placeholder function
}

void runTests(boost::ptr_vector<Test> &tests, std::string test_name)
{
	for( auto it = tests.begin(); it != tests.end(); ++it ){
		if(it->name == test_name)
		{
			std::cout << it->description << "..." << std::endl;
			it->run();
			if( it->excpt_occurred ){ 
				std::cout << "EXCEPTION" << std::endl;
			}else if( it->passed ){ 
				std::cout << "PASS" << std::endl;
			}
			else{ 
				std::cout << "FAIL" << std::endl;
			}
		}
	}
}
void runTests(boost::ptr_vector<Test> &tests){

	int num_passed = 0;
	int num_failed = 0;
	int num_run = 0;
	int num_exceptions = 0;

	boost::ptr_vector<Test>::iterator it = tests.begin();
	for( ; it != tests.end(); ++it ){
		std::cout << it->description << "..." << std::endl;
		it->run();
		if( it->excpt_occurred ){ 
			num_exceptions++;
			std::cout << "EXCEPTION" << std::endl;
		}else if( it->passed ){ 
			num_passed ++;
			std::cout << "PASS" << std::endl;
		}
		else{ 
			std::cout << "FAIL" << std::endl;
			num_failed++;
		}
		num_run++;
	}

	std::cout << std::endl << "----------------------------------------------" << std::endl;
	std::cout << num_run << " Tests Run" << std::endl;
	std::cout << num_passed << " Passed, ";
	std::cout << num_failed << " Failed, ";
	std::cout << num_exceptions << " Exceptions";
	std::cout << std::endl << "----------------------------------------------" << std::endl;

}

int main(int argc, char *argv[])
{
	boost::ptr_vector<Test> tests;
    MPI_Init( &argc, &argv );
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump == 1 ){
//		tests.push_back( new EMTestMiniBatchSelection() );
//		tests.push_back( new NNParamsTestBiasIndexes() );
//		tests.push_back( new NNParamsTestDropout());
		tests.push_back( new NNParamsTestComputeAndAccumulateGradient() );
//		tests.push_back( new NNParamsTestComputeUnweightedGradients() );
//		tests.push_back( new NNParamsTestComputeDeltas() );
//		tests.push_back( new NNParamsTestComputeTransitionThetas() );
//		tests.push_back( new NNParamsTestSaveAndLoadFromFile() );
//		//tests.push_back( new FragGenTestMaxRingBreaks() );
//		tests.push_back( new FVFragGraphSaveAndLoadState() );
//		tests.push_back( new FragGenTestDisallowDetourTransitions() );
//		tests.push_back( new MspReaderTest() );
//		tests.push_back( new ParamsTestComputeTransitionThetas() );
//		tests.push_back( new ParamsTestComputeAndAccumulateGradient());
//        tests.push_back( new MspReaderMultipleEnergiesTest() );
//        tests.push_back( new FragGenTestPositiveEIAlkane() );
//        //tests.push_back( new FragGenTestPositiveEINistExceptions() );
//        tests.push_back( new FragGenTestCasesFromGross() ); // NOTE THIS IS NOT WORKING
//        tests.push_back( new FragGenTestPositiveEIAndESIDegreeLpBonding() );
//        tests.push_back( new FragGenTestPositiveEIOxygenAromatic() );
        //tests.push_back( new FragGenTestPositiveEITriple() );
        //tests.push_back( new FragGenTestPositiveEISplitCharge() );
        //tests.push_back( new EMTestSingleEnergySelfProduction() );
//        tests.push_back( new EMTestNNSingleEnergySelfProduction() );
        //tests.push_back( new EMTestSelfProduction() );
        //tests.push_back( new EMTestSingleEnergyIsotopeSelfProduction() );
	}
	if( mpi_nump > 2 ){
		tests.push_back( new CommsTestSetMasterUsedIdxs() );
		tests.push_back( new CommsTestCollectQInMaster() );
		tests.push_back( new CommsTestCollectGradsInMaster() );
		tests.push_back( new CommsTestBroadcastParams() );
		tests.push_back( new EMTestMultiProcessor() );
	}

	// run one specific test
	if (argc == 2)
	{
		runTests( tests ,  argv[1]);
	}
	else
	{
		runTests( tests );
	}

	MPI_Finalize();

	return(0);    
}



