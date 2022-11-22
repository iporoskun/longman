#pragma once
#include <chrono>
#include <concepts>
#include "named_type.hpp"


namespace iporoskun::longman {

namespace detail {
struct latitude_tag;
struct longitude_tag;
struct height_tag;
} // namespace detail

template<std::floating_point FloatingType = double>
using latitude = detail::named_type<FloatingType, detail::latitude_tag>;

template<std::floating_point FloatingType = double>
using longitude = detail::named_type<FloatingType, detail::longitude_tag>;

template<std::floating_point FloatingType = double>
using height = detail::named_type<FloatingType, detail::height_tag>;

template<std::floating_point FloatingType = double>
struct position {
  latitude<FloatingType> latitude; /* deg */
  longitude<FloatingType> longitude; /* deg */
  height<FloatingType> height; /* msl_orthometric_height, meter and also cm*/
};

} // namespace iporoskun::longman