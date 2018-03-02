
macro(add_psc_test name)
  if(GTEST_FOUND)
    add_executable(${name} ${name}.cxx)
    target_link_libraries(${name} psc GTest::GTest GTest::Main)
    gtest_discover_tests(${name})
  endif()
endmacro()

add_psc_test(vec3_unittest)
add_psc_test(test_grid)
