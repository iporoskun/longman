# C++ implementation of Longman tidal equation

[![ci](https://github.com/iporoskun/longman/actions/workflows/ci.yml/badge.svg)](https://github.com/iporoskun/longman/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/iporoskun/longman/branch/main/graph/badge.svg)](https://codecov.io/gh/iporoskun/longman)
[![Language grade: C++](https://img.shields.io/lgtm/grade/cpp/github/cpp-best-practices/gui_starter_template)](https://lgtm.com/projects/g/cpp-best-practices/gui_starter_template/context:cpp)
[![CodeQL](https://github.com/iporoskun/longman/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/iporoskun/longman/actions/workflows/codeql-analysis.yml)

## Getting Started

The time point for the calculation must be given as UTC (+0h). To compute the acceleration at your local time you can do the following: 

```cpp
using namespace std::literals::chrono_literals;
const auto now = std::chrono::system_clock::now();
const auto now_as_utc = now - 1h; // -1 hour for UTC+1
// or directly use 
// const auto now_as_utc = std::chrono::utc_clock::now();
const auto acceleration = longman(my_position, now_as_utc);
```

## Useful references

Following sources can be referred for further details on Longman equations:

* [This publication](https://sbgf.org.br/revista/index.php/rbgf/article/viewFile/793/416) was used as a reference for mathematical details and parameters of the implementation.
* The source code is the C++ port of original MATLAB implementation by Olga Bjelotomic Orsulic and Matej Varga      
* Constants from [Advances in Geophysical Methods Applied to Forensic Investigations](https://shorturl.at/azJ49)


## Used libraries

* [Catch2](https://github.com/catchorg/Catch2/) for unit tests