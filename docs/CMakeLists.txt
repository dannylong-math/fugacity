set(FUGACITY_SPHINX_VENV "${PROJECT_SOURCE_DIR}/.dependencies/docs")

find_program(SPHINX_BUILD_EXECUTABLE
    NAMES sphinx-build
    HINTS "${FUGACITY_SPHINX_VENV}/bin"
    NO_DEFAULT_PATH
)

if(NOT SPHINX_BUILD_EXECUTABLE)
    find_program(SPHINX_BUILD_EXECUTABLE NAMES sphinx-build)
endif()

if(NOT SPHINX_BUILD_EXECUTABLE)
    message(FATAL_ERROR
        "sphinx-build was not found. Install the documentation dependencies with:\n"
        "  python3 -m venv .dependencies/docs\n"
        "  .dependencies/docs/bin/pip install -r docs/requirements.txt"
    )
endif()

file(GLOB_RECURSE FUGACITY_DOCUMENTATION_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.rst"
)
file(GLOB_RECURSE FUGACITY_PUBLIC_HEADERS CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
)

set(FUGACITY_DOCUMENTATION_OUTPUT "${CMAKE_BINARY_DIR}/html")
set(FUGACITY_DOCUMENTATION_DOCTREES "${CMAKE_BINARY_DIR}/doctrees")

add_custom_target(docs
    COMMAND
        "${SPHINX_BUILD_EXECUTABLE}"
        --builder html
        --doctree-dir "${FUGACITY_DOCUMENTATION_DOCTREES}"
        --fresh-env
        --fail-on-warning
        --keep-going
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${FUGACITY_DOCUMENTATION_OUTPUT}"
    DEPENDS ${FUGACITY_DOCUMENTATION_SOURCES} ${FUGACITY_PUBLIC_HEADERS}
    COMMENT "Building Fugacity documentation with Sphinx"
    VERBATIM
)

message(STATUS "Sphinx executable: ${SPHINX_BUILD_EXECUTABLE}")
message(STATUS "Documentation output: ${FUGACITY_DOCUMENTATION_OUTPUT}/index.html")
