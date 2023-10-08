# FindRDKit.cmake
# Placed in the public domain by NextMove Software in 2013
# Modified by Felicity Allen, August 2013
# Modified by Fei Wang,June 2017
# Try to find RDKit headers and libraries
# Defines:
#
#  RDKIT_FOUND - system has RDKit
#  RDKIT_INCLUDE_DIR - the RDKit include directory
#  RDKIT_INCLUDE_EXT_DIR - the RDKit external directory when including Inchi support
#  RDKIT_LIBRARIES - Link these to use RDKit

if (RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
    # in cache already or user-specified
    set(RDKIT_FOUND TRUE)

else ()

    if (NOT RDKIT_INCLUDE_DIR OR NOT RDKIT_INCLUDE_EXT_DIR)
        if (WIN32)
            find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
                    PATHS
                    ${RDKIT_DIR}\\Code
                    $ENV{RDKIT_INCLUDE_DIR}
                    $ENV{RDKIT_INCLUDE_PATH}
                    $ENV{RDKIT_BASE}\\Code
                    $ENV{RDBASE}\\Code
                    C:\\RDKit\\include
                    C:\\RDKit\\Code
                    )
            find_path(RDKIT_INCLUDE_EXT_DIR INCHI-API/inchi.h
                    PATHS
                    ${RDKIT_DIR}\\External
                    $ENV{RDKIT_INCLUDE_EXT_DIR}
                    $ENV{RDKIT_INCLUDE_EXT_PATH}
                    $ENV{RDKIT_BASE}\\External
                    $ENV{RDBASE}\\External
                    C:\\RDKit\\include
                    C:\\RDKit\\External
                    )
        else ()
            find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
                    PATHS
                    ${RDKIT_DIR}/Code
                    $ENV{RDKIT_INCLUDE_DIR}
                    $ENV{RDKIT_INCLUDE_PATH}
                    $ENV{RDKIT_BASE}/Code
                    $ENV{RDBASE}/Code
                    /usr/include/rdkit
                    /usr/local/include/rdkit
                    /usr/local/rdkit/include
                    /usr/local/rdkit/Code
                    /opt/rdkit/Code
                    /opt/rdkit/include
                    $ENV{RDBASE}/include/rdkit/External
                    ~/rdkit/Code
                    )
            find_path(RDKIT_INCLUDE_EXT_DIR INCHI-API/inchi.h
                    PATHS
                    ${RDKIT_DIR}/External
                    $ENV{RDKIT_INCLUDE_EXT_DIR}
                    $ENV{RDKIT_INCLUDE_EXT_PATH}
                    $ENV{RDKIT_BASE}/External
                    $ENV{RDBASE}/External
                    $ENV{RDBASE}/include/rdkit
                    /usr/include/rdkit
                    /usr/local/include/rdkit
                    /usr/local/rdkit/include
                    /usr/local/rdkit/Code
                    /opt/rdkit/Code
                    /opt/rdkit/include
                    /opt/rdkit/External
                    ~/rdkit/Code
                    /usr/include/rdkit/External
                    /usr/local/include/rdkit/External
                    /usr/local/rdkit/External
                    /usr/local/rdkit/External
                    ~/rdkit/External
                    )
        endif ()
        if (RDKIT_INCLUDE_DIR)
            message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_DIR}")
        else()
            message(STATUS "Can NOT Found RDKit include files at ${RDKIT_INCLUDE_DIR}")
        endif ()
        if (RDKIT_INCLUDE_EXT_DIR)
            message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_EXT_DIR}")
        else()
            message(STATUS "Can NOT Found RDKit include files at ${RDKIT_INCLUDE_EXT_DIR}")
        endif ()
    endif ()

    if (NOT RDKIT_LIBRARIES)
        if (NOT RDKIT_LIBRARY_DIR)
            find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers
                    PATHS
                    ${RDKIT_DIR}/lib
                    $ENV{RDKIT_LIB_DIR}
                    $ENV{RDKIT_LIB_PATH}
                    $ENV{RDKIT_LIBRARIES}
                    $ENV{RDKIT_BASE}/lib
                    $ENV{RDBASE}/lib
                    /usr/local/rdkit/lib
                    ~/rdkit/lib
                    /opt/rdkit/lib
                    $ENV{LD_LIBRARY_PATH}
                    )
        else ()
            find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers HINTS ${RDKIT_LIBRARY_DIR})
        endif ()
        if (FILEPARSERS_LIB)
            GET_FILENAME_COMPONENT(RDKIT_LIBRARY_DIR ${FILEPARSERS_LIB} PATH)
            message(STATUS "Found RDKit libraries at ${RDKIT_LIBRARY_DIR}")

            # Note that the order of the following libraries is significant!!
            find_library(SMILESPARSE_LIB NAMES SmilesParse RDKitSmilesParse
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(DEPICTOR_LIB NAMES Depictor RDKitDepictor
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(DESCIPTORS_LIB NAMES Descriptors RDKitDescriptors
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(CHEMTRANS_LIB NAMES ChemTransforms RDKitChemTransforms
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(GRAPHMOL_LIB NAMES GraphMol RDKitGraphMol
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(RDGEOMETRYLIB_LIB NAMES RDGeometryLib RDKitRDGeometryLib
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(RDGENERAL_LIB NAMES RDGeneral RDKitRDGeneral
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(SUBSTRUCT_LIB NAMES SubstructMatch RDKitSubstructMatch
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(GASTEIGER_LIB NAMES PartialCharges RDKitPartialCharges
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(DATASTRUCT_LIB NAMES DataStructs RDKitDataStructs
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(SUBGRAPH_LIB NAMES Subgraphs RDKitSubgraphs
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(FINGERPRINT_LIB NAMES Fingerprints RDKitFingerprints
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(INCHI_LIB NAMES Inchi RDKitInchi
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(RDINCHI_LIB NAMES RDInchiLib RDKitRDInchiLib
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(OPT NAMES Optimizer RDKitOptimizer
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(FF NAMES ForceField RDKitForceField
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(FFHELP NAMES ForceFieldHelpers RDKitForceFieldHelpers
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(CATALOG NAMES Catalogs RDKitCatalogs
                    HINTS ${RDKIT_LIBRARY_DIR})
            find_library(FRAGCAT NAMES FragCatalog RDKitFragCatalog
                    HINTS ${RDKIT_LIBRARY_DIR})


            set(RDKIT_LIBRARIES ${FILEPARSERS_LIB} ${SMILESPARSE_LIB}
                    ${DEPICTOR_LIB} ${CHEMTRANS_LIB} ${GRAPHMOL_LIB} ${RDGEOMETRYLIB_LIB}
                    ${RDGENERAL_LIB} ${SUBSTRUCT_LIB} ${GASTEIGER_LIB}
                    ${DATASTRUCT_LIB} ${SUBGRAPH_LIB} ${FINGERPRINT_LIB}
                    ${INCHI_LIB} ${RDINCHI_LIB} ${DESCIPTORS_LIB} ${OPT} ${FF} ${FFHELP} ${CATALOG} ${FRAGCAT})
        endif ()
        if (RDKIT_LIBRARIES)
            message(STATUS "Found RDKit library files at ${RDKIT_LIBRARIES}")
        endif ()
    endif ()

    if (RDKIT_INCLUDE_DIR AND RDKIT_INCLUDE_EXT_DIR AND RDKIT_LIBRARIES)
        set(RDKIT_FOUND TRUE)
    endif ()

    mark_as_advanced(RDINCHI_LIB INCHI_LIB GASTEIGER_LIB SUBSTRUCT_LIB RDGENERAL_LIB RDGEOMETRYLIB_LIB GRAPHMOL_LIB DEPICTOR_LIB SMILESPARSE_LIB FILEPARSERS_LIB)
    mark_as_advanced(RDKIT_INCLUDE_DIR RDKIT_INCLUDE_EXT_DIR RDKIT_LIBRARIES)
endif ()
