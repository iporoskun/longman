cmake_minimum_required(VERSION 3.15...3.23)

project(CmakeConfigPackageTests LANGUAGES CXX)

# ---- Test as standalone project the exported config package ----
if(PROJECT_IS_TOP_LEVEL OR TEST_INSTALLED_VERSION)
  enable_testing()

  find_package(longman CONFIG REQUIRED) # for tests, project_options, ...

  if(NOT TARGET longman::project_options)
    message(FATAL_ERROR "Requiered config package not found!")
    return() # be strictly paranoid for Template Janitor github action! CK
  endif()
endif()

# ---- Dependencies ----
find_package(Catch2 CONFIG REQUIRED)

include(Catch)

add_library(catch_main OBJECT catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch2)
target_link_libraries(catch_main PRIVATE longman::project_options)

set_property(TARGET catch_main PROPERTY FOLDER "test")

function(IPOROSKUN_LONGMAN_add_test path name)
  add_executable(${name} ${path})
  target_link_libraries(${name} PUBLIC
    longman::longman
    longman::project_options
    longman::project_warnings
    catch_main
  )
  catch_discover_tests(${name}
    TEST_PREFIX "${name}.."
    EXTRA_ARGS
    -s
    --reporter=xml
    --out=${name}.xml)

  target_compile_definitions(${name} PRIVATE
  )

  set_property(TARGET ${name} PROPERTY FOLDER "test")
endfunction()

file(GLOB_RECURSE IPOROSKUN_LONGMAN_test_list "${CMAKE_CURRENT_LIST_DIR}/*.cpp")
LIST(REMOVE_DUPLICATES IPOROSKUN_LONGMAN_test_list)
LIST(REMOVE_ITEM IPOROSKUN_LONGMAN_test_list "${CMAKE_CURRENT_LIST_DIR}/catch_main.cpp")

FOREACH(file_path ${IPOROSKUN_LONGMAN_test_list})
  get_filename_component(test_name ${file_path} NAME_WE)
  IPOROSKUN_LONGMAN_add_test(${file_path} ${test_name})
ENDFOREACH()
