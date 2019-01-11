##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-code/CMakeLists.txt
#
# Author: Fei Wang
# Created: Jan 2019
##########################################################################
set ( SRC_FILES main.cpp)

add_executable ( cfm-boost-test ${SRC_FILES} )
target_link_libraries ( cfm-boost-test cfm-code ${REQUIRED_LIBS} )

# declares a test with our executable
#add_test(NAME test1 COMMAND test_executable)