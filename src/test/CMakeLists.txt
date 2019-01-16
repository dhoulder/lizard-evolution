include_directories(..)

# FIXME doesn't have automatic dependency on ./dreadds
add_test(NAME basic
         COMMAND dreadds-bin -c 1/test.conf -o 1/out)

add_executable(test-speciation test-speciation.cpp )
target_link_libraries(test-speciation  dreadds)
set_property(TARGET test-speciation
             PROPERTY CXX_STANDARD 11)
add_test(NAME speciation
         COMMAND test-speciation)

add_test(NAME Setup COMMAND ./setup.sh)
add_test(NAME Cleanup COMMAND ./cleanup.sh)

set_tests_properties(Setup PROPERTIES FIXTURES_SETUP BASIC)
set_tests_properties(Cleanup PROPERTIES FIXTURES_CLEANUP BASIC)
set_tests_properties(basic  PROPERTIES FIXTURES_REQUIRED BASIC)
