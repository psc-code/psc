
add_subdirectory(libpsc)

add_executable(psc_bubble_yz psc_bubble_yz.cxx)
target_include_directories(psc_bubble_yz PRIVATE include)
target_link_libraries(psc_bubble_yz psc)

#target_compile_features(mylib PUBLIC cxx_std_11)

