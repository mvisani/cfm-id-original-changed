##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-code/CMakeLists.txt
#
# Author: Felicity Allen, Fei Wang
# Created: August 2013
# Last Modified: Sep 2017
##########################################################################

set(BASE_HEADERS FragmentGraph.h Feature.h FeatureVector.h FeatureCalculator.h Config.h Spectrum.h
        FragmentGraphGenerator.h FragmentTreeNode.h Inference.h MILP.h Isotope.h
        MolData.h MspReader.h Param.h NNParam.h Util.h Message.h Identifier.h Comparators.h
        FunctionalGroups.h)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Features FEATURES_SRC_DIR)

set(TRAIN_HEADERS Comms.h ModelBase.h EmModel.h EmNNModel.h Solver.h)

set(BASE_SRC_FILES FragmentGraph.cpp
        Feature.cpp
        FeatureVector.cpp
        FeatureCalculator.cpp
        ${FEATURES_SRC_DIR}
        Config.cpp
        FragmentGraphGenerator.cpp
        FragmentTreeNode.cpp
        FunctionalGroups.cpp
        Identifier.cpp
        Inference.cpp
        Isotope.cpp
        Message.cpp
        MILP.cpp
        MolData.cpp
        MspReader.cpp
        Param.cpp
        NNParam.cpp
        Spectrum.cpp
        Util.cpp
        Comparators.cpp)

set(TRAIN_SRC_FILES Comms.cpp
        ModelBase.cpp
        EmModel.cpp
        EmNNModel.cpp
        Solver.cpp)

if (INCLUDE_TRAIN)
    set(SRC_FILES ${BASE_SRC_FILES} ${TRAIN_SRC_FILES})
    set(HEADER_FILES ${BASE_HEADERS} ${TRAIN_HEADERS})
    set(REQUIRED_LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES}
            ${LPSOLVE_LIBRARY} ${MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
else ()
    set(SRC_FILES ${BASE_SRC_FILES})
    set(HEADER_FILES ${BASE_HEADERS})
    set(REQUIRED_LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES}
            ${LPSOLVE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
endif ()

add_library(cfm-code STATIC ${SRC_FILES} ${HEADER_FILES})
target_link_libraries(cfm-code ${REQUIRED_LIBS})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ISOTOPE.DAT DESTINATION ${CMAKE_BINARY_DIR}/bin)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ISOTOPE.DAT
        DESTINATION ${CFM_OUTPUT_DIR})