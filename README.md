# C++ implementation of Longman tidal equation

[![ci](https://github.com/iporoskun/longman/actions/workflows/ci.yml/badge.svg)](https://github.com/iporoskun/longman/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/iporoskun/longman/branch/main/graph/badge.svg)](https://codecov.io/gh/iporoskun/longman)
[![Language grade: C++](https://img.shields.io/lgtm/grade/cpp/github/cpp-best-practices/gui_starter_template)](https://lgtm.com/projects/g/cpp-best-practices/gui_starter_template/context:cpp)
[![CodeQL](https://github.com/iporoskun/longman/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/iporoskun/longman/actions/workflows/codeql-analysis.yml)


## Getting started

### How to use it

`ip::longman` is a header only library. You can simply take the file `include\iporoskun\longman.hpp` and include it into your project.


### Usage example

For convenience, we declare namespace alias `ip`
```
namespace ip = iporoskun::longman;
```

First, we need to specify the geographical position of the point on earth for which the computation should be performed. The position is descripted by latitude (in degrees), longitude (in degrees), and the height above mean sea level (in meters). 

As an example, we will use the location of Brandenburger Tor (Berlin): N 52.51628 E 13.377702 and 38 m height. 

```cpp
const auto position = ip::position {
    ip::latitude_t<>(52.51628),
    ip::longitude_t<>(13.377702),
    ip::height_t<>(38.)
};
```
To specify the position, the library provides the class template position as well as three class templates for the position components. 

Second, we need a time point for which we want to perform the computation. The time point for the calculation must be given as UTC time. Here I will just use the current time on my machine:
```cpp
const auto current_time = std::chrono::utc_clock::now();
```

The function template `ip::longman` takes the position and time point, and computes the tidal acceleration in meters per second squared. The results is returned as a floating point number: 
```cpp
const double delta_acceleration = ip::longman(position, current_time);
```

### Details

The function template `ip::longman` has the following signature:
```cpp
template<std::floating_point FloatingType, class DurationType, class TimePoint>
[[nodiscard]] inline auto ip::longman(const ip::position<FloatingType>& position, const TimePoint& utc_time)
	-> FloatingType
```

As the first argument, the used floating point type can be specified; any of ~~`float`~~, `double`, and `long double`. 

The time point (`TimePoint`) can be provided by (theoretically) *any* clock type; although it was only tested for `std::chrono::system_clock` and `std::chrono::gps_clock`. 
> No matter what clock is used, the time must always be in UTC format. 


The duration (`DurationType`) can be anything from `std::chrono::nanoseconds` to `std::chrono::minutes`. (Usage of ` std::chrono::hours` gave slightly different results, probably due to rounding of the minutes to whole hours.)


## Benchmark 

Benchmark results from running in release mode:
```
Run on (8 X 1500.44 MHz CPU s)
CPU Caches:
  L1 Data 48 KiB (x4)
  L1 Instruction 32 KiB (x4)
  L2 Unified 512 KiB (x4)
  L3 Unified 6144 KiB (x1)
--------------------------------------------------------
Benchmark              Time             CPU   Iterations
--------------------------------------------------------
longman_bench       1122 ns         1116 ns       560000
```


## Useful references

Following sources can be referred for further details on Longman equations:

* [This publication](https://sbgf.org.br/revista/index.php/rbgf/article/viewFile/793/416) was used as a reference for mathematical details and parameters of the implementation.
* The source code is the C++ port of original MATLAB implementation by Olga Bjelotomic Orsulic and Matej Varga      
* Constants from [Advances in Geophysical Methods Applied to Forensic Investigations](https://shorturl.at/azJ49)


## Used libraries

* [Catch2](https://github.com/catchorg/Catch2/) for unit tests
* [google::benchmark](https://github.com/google/benchmark) for benchmarking