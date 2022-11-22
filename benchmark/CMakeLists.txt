
add_executable(longman_bench longman_bench.cpp)
target_link_libraries(longman_bench longman::longman benchmark::benchmark_main)
