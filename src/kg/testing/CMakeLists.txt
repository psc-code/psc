
macro(add_kg_test name src)
  add_executable(${name})
  target_gtensor_sources(${name} PRIVATE ${src})
  target_link_libraries(${name} GTest::gtest GTest::gtest_main kg psc)
  if (USE_GTEST_DISCOVER_TESTS)
    gtest_discover_tests(${name})
  else()
    gtest_add_tests(TARGET ${name})
  endif()
endmacro()

add_kg_test(TestVec3 TestVec3.cxx)
add_kg_test(TestSArray TestSArray.cxx)

if (PSC_HAVE_ADIOS2)
  add_kg_test(TestIOAdios2 io/TestIOAdios2.cxx)
  add_kg_test(TestIO io/TestIO.cxx)
endif()
