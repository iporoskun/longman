#pragma once
#include <chrono>
#include "named_type.hpp"

#ifndef IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE
#  define IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE double
#endif

namespace iporoskun::longman {

using floating_t = IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE;
using duration_t = std::chrono::seconds;

using latitude = details::named_type<floating_t, struct latitude_tag>;
using longitude = details::named_type<floating_t, struct longitude_tag>;
using height = details::named_type<floating_t, struct height_tag>;

struct position {
  latitude latitude; /* deg */
  longitude longitude; /* deg */
  height height; /* msl_orthometric_height, meter and also cm*/
};

} // namespace iporoskun::longman