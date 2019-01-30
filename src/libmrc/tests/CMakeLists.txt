
link_libraries(mrc)

add_executable(test_mrc_obj test_mrc_obj.c)
target_compile_features(test_mrc_obj
  PRIVATE
    c_std_99
)

add_executable(test_io test_io.c)
target_compile_features(test_io
  PRIVATE
    c_std_99
)

