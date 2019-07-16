
set(INSTALL_GTEST OFF CACHE INTERNAL "")
set(BUILD_GMOCK OFF CACHE INTERNAL "")

add_subdirectory(googletest EXCLUDE_FROM_ALL)

add_library(GTest::GTest ALIAS gtest)
add_library(GTest::Main ALIAS gtest_main)
