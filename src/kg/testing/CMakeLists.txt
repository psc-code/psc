
macro(add_kg_test name src)
  add_executable(${name} ${src})
  target_link_libraries(${name} GTest::GTest GTest::Main kg)
  gtest_discover_tests(${name})
endmacro()

add_kg_test(TestVec3 TestVec3.cxx)
  
if (PSC_HAVE_ADIOS2)
  add_kg_test(TestIOAdios2 io/TestIOAdios2.cxx)
  add_kg_test(TestIO io/TestIO.cxx)
endif()
