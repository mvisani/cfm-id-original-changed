##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-code/CMakeLists.txt
#
# Author: Fei Wang
# Created: Jan 2019
##########################################################################
set(SRC_FILES main.cpp features_tests/FeaturesTestsUtil.cpp features_tests/FeaturesTestRootAtom.cpp features_tests/FeaturesTestsUtil.h features_tests/FeaturesTestsInitClass.cpp features_tests/FeaturesTestBreakAtomPair.cpp features_tests/FeaturesTestRootPairs.cpp features_tests/FeaturesTestRootTriples.cpp features_tests/FeaturesTestGasteigerCharges.cpp features_tests/FeaturesTestHydrogenMovement.cpp features_tests/FeaturesTestFunctionalGroups.cpp features_tests/FeaturesTestExtraFunctionalGroups.cpp features_tests/FeaturesTestFunctionalGroupsRootOnly.cpp features_tests/FeaturesTestRadicalFeatures.cpp features_tests/FeaturesTestExtraRingFeatures.cpp features_tests/FeaturesTestRootMMFFAtomType.cpp features_tests/FeaturesTestNeighbourMMFFAtomType.cpp features_tests/FeaturesTestBrokenOrigBondType.cpp features_tests/FeaturesTestQuadraticFeatures.cpp features_tests/FeaturesTestLength.cpp features_tests/FeaturesTestMatrixFP.cpp features_tests/FeaturesTestIonicFeatures.cpp)

add_executable ( cfm-boost-test ${SRC_FILES} )
target_link_libraries ( cfm-boost-test cfm-code ${REQUIRED_LIBS})

# declares a test with our executable
#enable_testing()
install ( TARGETS cfm-boost-test
        DESTINATION ${CFM_OUTPUT_DIR} )

install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_data
        DESTINATION ${CFM_OUTPUT_DIR}/tests )