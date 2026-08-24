find_package(Doxygen 1.9.8 REQUIRED)
find_program(FUGACITY_NODE_EXECUTABLE NAMES node REQUIRED)
find_program(FUGACITY_NPM_EXECUTABLE NAMES npm REQUIRED)

execute_process(
    COMMAND "${FUGACITY_NODE_EXECUTABLE}" --version
    OUTPUT_VARIABLE FUGACITY_NODE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)
string(REGEX REPLACE "^v" "" FUGACITY_NODE_VERSION "${FUGACITY_NODE_VERSION}")
if(FUGACITY_NODE_VERSION VERSION_LESS 22.12)
    message(FATAL_ERROR
        "Sourcey requires Node.js 22.12 or newer; found ${FUGACITY_NODE_VERSION}"
    )
endif()

file(GLOB FUGACITY_DOCUMENTATION_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.md"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.dox"
)
file(GLOB_RECURSE FUGACITY_PUBLIC_HEADERS CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
)

set(FUGACITY_DOXYGEN_XML "${PROJECT_SOURCE_DIR}/build/doxygen/xml/index.xml")
set(FUGACITY_SOURCEY_OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/dist/index.html")

add_custom_target(docs
    COMMAND "${CMAKE_COMMAND}" -E make_directory
            "${PROJECT_SOURCE_DIR}/build/doxygen"
    COMMAND "${DOXYGEN_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/Doxyfile"
    COMMAND "${FUGACITY_NPM_EXECUTABLE}" ci --prefix
            "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMAND "${FUGACITY_NPM_EXECUTABLE}" run --prefix
            "${CMAKE_CURRENT_SOURCE_DIR}" build
    BYPRODUCTS
        "${FUGACITY_DOXYGEN_XML}"
        "${FUGACITY_SOURCEY_OUTPUT}"
    DEPENDS
        "${PROJECT_SOURCE_DIR}/Doxyfile"
        "${CMAKE_CURRENT_SOURCE_DIR}/package.json"
        "${CMAKE_CURRENT_SOURCE_DIR}/package-lock.json"
        "${CMAKE_CURRENT_SOURCE_DIR}/sourcey.config.ts"
        "${CMAKE_CURRENT_SOURCE_DIR}/filter-doxygen-xml.mjs"
        "${CMAKE_CURRENT_SOURCE_DIR}/render-math.mjs"
        ${FUGACITY_DOCUMENTATION_SOURCES}
        ${FUGACITY_PUBLIC_HEADERS}
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    COMMENT "Building Fugacity documentation with Doxygen and Sourcey"
    VERBATIM
)

message(STATUS "Doxygen executable: ${DOXYGEN_EXECUTABLE}")
message(STATUS "Node.js executable: ${FUGACITY_NODE_EXECUTABLE} (${FUGACITY_NODE_VERSION})")
message(STATUS "Documentation output: ${FUGACITY_SOURCEY_OUTPUT}")
