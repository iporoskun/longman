#pragma once
#include <chrono>
#include <concepts>
#include "named_type.hpp"

#ifndef IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE
#  define IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE double
#endif

// if needed, define to disable: IPOROSKUN_LONGMAN_DISABLE_RUNTIME_CHECKS

namespace iporoskun::longman {

namespace detail {
using floating_t = IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE;

struct latitude_tag;
struct longitude_tag;
struct height_tag;
} // namespace detail

template<std::floating_point FloatingType = detail::floating_t>
using latitude = detail::named_type<FloatingType, detail::latitude_tag>;

template<std::floating_point FloatingType = detail::floating_t>
using longitude = detail::named_type<FloatingType, detail::longitude_tag>;

template<std::floating_point FloatingType = detail::floating_t>
using height = detail::named_type<FloatingType, detail::height_tag>;

template<std::floating_point FloatingType = detail::floating_t>
struct position {
  latitude<FloatingType> latitude; /* deg */
  longitude<FloatingType> longitude; /* deg */
  height<FloatingType> height; /* msl_orthometric_height, meter and also cm*/
};

} // namespace iporoskun::longman