# List of sources
list(APPEND sources_list ddcosmo.f90 sphere_lebedev_rule.f90)

add_library(ddpcm OBJECT ${sources_list})
set_target_properties(ddpcm PROPERTIES POSITION_INDEPENDENT_CODE 1)
