
find_package(fmt CONFIG)

# Generic test that uses conan libs
add_executable(intro main.cpp)
target_link_libraries(
  intro
  PUBLIC project_options project_warnings
  PRIVATE fmt::fmt)

target_include_directories(intro PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include")
