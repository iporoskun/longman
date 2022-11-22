#include <iporoskun/longman.hpp>
#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <format>

using namespace std::chrono_literals;
namespace ip = iporoskun::longman;
using floating_t = long double;

const auto folder_path = std::filesystem::current_path().string();

const auto position = ip::position{ ip::latitude_t<floating_t>(52.51628),
									ip::longitude_t<floating_t>(13.377702),
									ip::height_t<floating_t>(38.) };
const auto start_time =
  std::chrono::sys_days{ 2019y / 10 / 1d } + 0h + 0min + 0s;

const auto time_increment = 1min;
const auto duration = 72h;


int main() {

  std::stringstream stream;

  auto time = start_time;
  while (time - start_time < duration) {
	const auto accel = ip::longman<floating_t>(position, time);
	const auto formated_time = std::chrono::round<std::chrono::minutes>(time);
	stream << formated_time << "," << accel * 1e9 << '\n'; // << " nm/s^2\n";

	time += time_increment;
  }

  std::cout << stream.str();

  std::ofstream myfile;
  myfile.open("longman_result.csv");
  myfile << stream.rdbuf();
}