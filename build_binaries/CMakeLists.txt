CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

project(ilastik)

#############################################
#
#     install 'buildem' build scripts
#
#############################################

set (BUILDEM_DIR "None" CACHE TYPE STRING)

if (${BUILDEM_DIR} STREQUAL "None")
    message (FATAL_ERROR "ERROR: FlyEM build directory (for all downloads & builds) should be specified via -D BUILDEM_DIR=<path> on cmake command line.")
endif ()

get_filename_component(BUILDEM_REPO_DIR ../buildem ABSOLUTE)
if (NOT EXISTS ${BUILDEM_REPO_DIR}/python.cmake)
    message ("Downloading buildem repository ...")
    execute_process(COMMAND git clone https://github.com/martinsch/buildem.git ${BUILDEM_REPO_DIR}
                    WORKING_DIRECTORY ${BUILDEM_DIR}
                    RESULT_VARIABLE BUILDEM_INSTALL_FAILED)
    if(BUILDEM_INSTALL_FAILED)
        MESSAGE(FATAL_ERROR "Unable to clone buildem scripts from\n     'https://github.com/martinsch/buildem.git'\n(either git is not installed properly on your system, or the repository is not accessible). Please place buildem at\n    '${BUILDEM_REPO_DIR}'.")
    endif()
else()
    message ("Updating buildem repoository in ${BUILDEM_REPO_DIR}")
    execute_process(COMMAND git pull
                    WORKING_DIRECTORY ${BUILDEM_REPO_DIR}
                    RESULT_VARIABLE BUILDEM_UPDATE_FAILED)
    if(BUILDEM_UPDATE_FAILED)
        MESSAGE(FATAL_ERROR "Unable to update buildem repository at ${BUILDEM_REPO_DIR}")
    endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${BUILDEM_REPO_DIR})

#############################################
#
#     compile ilastik and its dependencies
#
#############################################

set(CPLEX_ROOT_DIR "" CACHE PATH "CPLEX root directory")

include(ilastik-gui)

set (PACKAGE_WITH_CPLEX OFF CACHE BOOL "Include CPLEX shared libraries in the final binary pacakge." FORCE)
if (PACKAGE_WITH_CPLEX)
	if (CPLEX_ROOT_DIR)
		install( TARGETS cplex-shared concert-shared ilocplex-shared
				 LIBRARY DESTINATION ${BUILDEM_DIR}/lib )
	else()
		message(FATAL_ERROR "Cannot package with CPLEX; CPLEX_ROOT_DIR is not valid")
	endif()
endif()

#############################################
#
#     install
#
#############################################

configure_file(${PROJECT_SOURCE_DIR}/relative-rpath.sh.in ${BUILDEM_DIR}/tmp/relative-rpath.sh @ONLY)
configure_file(${PROJECT_SOURCE_DIR}/list-dependencies.sh ${BUILDEM_DIR}/bin/list-dependencies.sh COPYONLY)

# compile the python sources so that we can potentially install into write-protected directories
install(CODE "execute_process(COMMAND ${PYTHON_EXE} ${PROJECT_SOURCE_DIR}/compile_all.py ${ilastik_SRC_DIR})")

# set relative RPATH in libraries and executable, or delete RPATH if necessary
install(CODE "execute_process(COMMAND ${BUILDEM_DIR}/tmp/relative-rpath.sh ${BUILDEM_DIR}/lib ${BUILDEM_DIR}/lib)")
install(CODE "execute_process(COMMAND ${BUILDEM_DIR}/tmp/relative-rpath.sh ${BUILDEM_DIR}/lib ${BUILDEM_DIR}/plugins)")
install(CODE "execute_process(COMMAND ${BUILDEM_DIR}/tmp/relative-rpath.sh delete ${ilastik_SRC_DIR})")
install(CODE "execute_process(COMMAND chrpath -r \\$ORIGIN/../lib ${BUILDEM_DIR}/bin/python2.7)")

INSTALL(PROGRAMS
    ${BUILDEM_DIR}/bin/ilastik_clusterized
    ${BUILDEM_DIR}/bin/ilastik_gui
    ${BUILDEM_DIR}/bin/ilastik_gui_test
    ${BUILDEM_DIR}/bin/ilastik_headless
    ${BUILDEM_DIR}/bin/ilastik_headless_test
    ${BUILDEM_DIR}/bin/setenv_ilastik_gui.sh
    ${BUILDEM_DIR}/bin/list-dependencies.sh
    ${BUILDEM_DIR}/bin/python   
    ${BUILDEM_DIR}/bin/python2   
    ${BUILDEM_DIR}/bin/python2.7    
  DESTINATION bin)

INSTALL(PROGRAMS ${PROJECT_SOURCE_DIR}/run_ilastik.sh
        DESTINATION bin/..)

INSTALL(FILES ${PROJECT_SOURCE_DIR}/qt.conf
        DESTINATION bin)

INSTALL(DIRECTORY ${ilastik_SRC_DIR}/ DESTINATION src/${ilastik_NAME} USE_SOURCE_PERMISSIONS)
INSTALL(DIRECTORY ${BUILDEM_DIR}/lib/ DESTINATION lib USE_SOURCE_PERMISSIONS)
INSTALL(DIRECTORY ${BUILDEM_DIR}/include/python2.7/ DESTINATION include/python2.7 USE_SOURCE_PERMISSIONS)
INSTALL(DIRECTORY ${BUILDEM_DIR}/plugins/ DESTINATION plugins USE_SOURCE_PERMISSIONS)

#############################################
#
#     create binary package
#
#############################################

SET(CPACK_PACKAGE_VENDOR "The ilastik Team")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "Interactive image labeling and segmentation toolkit."
)

# package version setup
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "beta")

SET(CPACK_PACKAGE_INSTALL_DIRECTORY     "${PROJECT_NAME}-1.0.beta")
SET(CPACK_STRIP_FILES TRUE)
SET(CPACK_PACKAGE_CONTACT "ilastik Team <team@ilastik.org>")
SET(CPACK_GENERATOR TGZ)

INCLUDE (CPack)
