#include <chrono>
#include <benchmark/benchmark.h>
#include <iporoskun/longman.hpp>

using namespace std::chrono;
namespace ip = iporoskun::longman;


const auto my_day = sys_days{ (year_month_day{ 2010y / 10 / 31d }) };
const auto my_time = my_day + 6h + 10min + 0s;

const auto position =
  ip::position<>{ ip::latitude_t<>(-22.733), ip::longitude_t<>(-90.50),
				  ip::height_t<>(0) };


static void longman_bench(benchmark::State& state) {
  for (auto _ : state) {
	[[maybe_unused]] const auto discard = ip::longman(position, my_time);
  }
}

BENCHMARK(longman_bench);
BENCHMARK_MAIN();