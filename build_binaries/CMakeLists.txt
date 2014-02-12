CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project (Foo)

include (ExternalProject)

############################################################################
# Check if BUILDEM_DIR has already been assigned.  If not, create a default.
set (BUILDEM_DIR "None" CACHE TYPE STRING)

message("BUILDEM_DIR = ${BUILDEM_DIR}")

if (${BUILDEM_DIR} STREQUAL "None")
    message (FATAL_ERROR "ERROR: BuildEM build directory (for all downloads & builds) should be specified via -D BUILDEM_DIR=<path> on cmake command line.")
endif ()

message ("BuildEM downloads and builds will be placed here: ${BUILDEM_DIR}")

############################################################################

############################################################################
# Download and install buildem, if it isn't already in BUILDEM_DIR.
set (BUILDEM_REPO_DIR ${BUILDEM_DIR}/src/buildem)
if (NOT EXISTS ${BUILDEM_REPO_DIR}/python.cmake)
    message ("Installing buildem repo...")
    ExternalProject_Add(buildem
        PREFIX              ${BUILDEM_DIR}
        GIT_REPOSITORY      https://github.com/martinsch/buildem.git
        #GIT_TAG            python3  # Example of tagged branch (see doc)
        UPDATE_COMMAND      ""
        PATCH_COMMAND       ""
        CONFIGURE_COMMAND   "" 
        BUILD_COMMAND       ""
        BUILD_IN_SOURCE     1
        INSTALL_COMMAND     ""
    )
    message ("\n**********************************************************\n")
    message ("\nAfter running make, you must re-run the cmake command once")
    message ("buildem has been downloaded!\n")
    message ("\n***********************************************************\n")
else ()
    ############################################################################

    # Use modules from the downloaded buildem
    set (CMAKE_MODULE_PATH ${BUILDEM_REPO_DIR})
    message("Using cmake modules from ${BUILDEM_REPO_DIR}")

    # Download, compile, test, and install dependencies
    # Note the auto-generated APP_DEPENDENCIES variable holds all included targets.
    include (ilastik-gui)

    add_custom_target (AppDependencies ALL
        DEPENDS ${APP_DEPENDENCIES}
        COMMENT "Installed all required software using buildem.")

    # Install Foo -- actual build commands should go below
    # add_executable(...)

############################################################################
endif()
