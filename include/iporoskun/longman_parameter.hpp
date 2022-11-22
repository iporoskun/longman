#pragma once
#include <chrono>
#include "named_type.hpp"

#ifndef IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE
#  define IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE double
#endif

// if needed, define to disable: IPOROSKUN_LONGMAN_DISABLE_RUNTIME_CHECKS

namespace iporoskun::longman {

using floating_t = IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE;
using duration_t = std::chrono::seconds;

namespace detail {
struct latitude_tag;
struct longitude_tag;
struct height_tag;
} // namespace detail

using latitude = detail::named_type<floating_t, detail::latitude_tag>;
using longitude = detail::named_type<floating_t, detail::longitude_tag>;
using height = detail::named_type<floating_t, detail::height_tag>;

struct position {
  latitude latitude; /* deg */
  longitude longitude; /* deg */
  height height; /* msl_orthometric_height, meter and also cm*/
};

} // namespace iporoskun::longman