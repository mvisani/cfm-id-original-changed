##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-feature-calculator/CMakeLists.txt
#
# Author: Fei Wang
# Created: August 2017
##########################################################################
 
set ( SRC_FILES  main.cpp )
                                           
add_executable ( cfm-feature-calculator ${SRC_FILES} )
target_link_libraries ( cfm-feature-calculator cfm-code ${REQUIRED_LIBS} )

install ( TARGETS cfm-feature-calculator 
          DESTINATION ${CFM_OUTPUT_DIR} )
          
          
          