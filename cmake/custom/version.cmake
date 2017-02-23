file(READ ${PROJECT_SOURCE_DIR}/README.md _readme)
string(REGEX MATCH "[0-9]\\.[0-9]\\.[0-9]" _version_string ${_readme})
string(REPLACE "." ";" _version_list ${_version_string})
list(GET _version_list 0 PROJECT_VERSION_MAJOR)
list(GET _version_list 1 PROJECT_VERSION_MINOR)
list(GET _version_list 2 PROJECT_VERSION_PATCH)

set(${PROJECT_NAME}_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
message(STATUS "${BoldGreen}PCMSolver v${${PROJECT_NAME}_VERSION}${ColourReset}")
