
macro(add_psc_test name)
  add_executable(${name} ${name}.cxx)
  target_link_libraries(${name} psc GTest::GTest GTest::Main)
  gtest_discover_tests(${name})
endmacro()

add_psc_test(test_grid)
