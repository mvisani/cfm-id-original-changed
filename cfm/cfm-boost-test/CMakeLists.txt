##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-code/CMakeLists.txt
#
# Author: Fei Wang
# Created: Jan 2019
##########################################################################
set(SRC_FILES main.cpp)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/features_tests FEATURES_TESTS_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/inference_tests INFERENCE_TESTS_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/fraggen_tests FRAGGEN_TESTS_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/util_tests UTIL_TESTS_SRC)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/learner_tests LEARNER_TESTS_SRC)

add_executable ( cfm-boost-test ${SRC_FILES} ${FEATURES_TESTS_SRC} ${INFERENCE_TESTS_SRC}
        ${FRAGGEN_TESTS_SRC} ${UTIL_TESTS_SRC} ${LEARNER_TESTS_SRC})
target_link_libraries ( cfm-boost-test cfm-code ${REQUIRED_LIBS})

# declares a test with our executable
#enable_testing()
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_data DESTINATION ${CMAKE_BINARY_DIR}/bin)

install ( TARGETS cfm-boost-test
        DESTINATION ${CFM_OUTPUT_DIR} )

install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_data
        DESTINATION ${CFM_OUTPUT_DIR}/test_data)