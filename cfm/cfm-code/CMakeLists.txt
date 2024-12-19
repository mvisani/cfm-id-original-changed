aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Features FEATURES_SRC_DIR)
set(BASE_HEADERS
    Comparators.h
    Config.h
    EmModel.h
    EmNNModel.h
    Feature.h
    FeatureCalculator.h
    FeatureVector.h
    FragmentGraph.h
    FragmentGraphGenerator.h
    FragmentTreeNode.h
    FunctionalGroups.h
    Identifier.h
    Inference.h
    Isotope.h
    MILP.h
    Message.h
    ModelBase.h
    MolData.h
    MspReader.h
    NNParam.h
    Param.h
    Solver.h
    Spectrum.h
    Util.h
    Version.h
)
set(BASE_SOURCES
    ${FEATURES_SRC_DIR}
    Comparators.cpp
    Config.cpp
    EmModel.cpp
    EmNNModel.cpp
    Feature.cpp
    FeatureCalculator.cpp
    FeatureVector.cpp
    FragmentGraph.cpp
    FragmentGraphGenerator.cpp
    FragmentTreeNode.cpp
    FunctionalGroups.cpp
    Identifier.cpp
    Inference.cpp
    Isotope.cpp
    MILP.cpp
    Message.cpp
    ModelBase.cpp
    MolData.cpp
    MspReader.cpp
    NNParam.cpp
    Param.cpp
    Solver.cpp
    Spectrum.cpp
    Util.cpp
)

add_library(cfm-code STATIC ${BASE_HEADERS} ${BASE_SOURCES})  # Include both headers and sources
set(REQUIRED_LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES} ${LPSOLVE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${OpenMP_CXX_FLAGS})
target_link_libraries(cfm-code ${REQUIRED_LIBS})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ISOTOPE.DAT DESTINATION ${CMAKE_BINARY_DIR}/bin)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ISOTOPE.DAT DESTINATION ${CFM_OUTPUT_DIR})
