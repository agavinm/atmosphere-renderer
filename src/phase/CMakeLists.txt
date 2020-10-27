set(MTS_PLUGIN_PREFIX "phasefunctions")

add_plugin(cha cha.cpp)
add_plugin(hg hg.cpp)
add_plugin(isotropic isotropic.cpp)
add_plugin(rayleigh rayleigh.cpp)

# Register the test directory
add_tests(${CMAKE_CURRENT_SOURCE_DIR}/tests)
