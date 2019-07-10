
macro(add_kg_test name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} GTest::GTest kg)
  gtest_discover_tests(${name})
endmacro()
  
if (PSC_HAVE_ADIOS2)
  add_kg_test(TestIOAdios2 io/TestIOAdios2.cxx)
  add_kg_test(TestIO io/TestIO.cxx)
endif()
