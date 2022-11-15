#pragma once
#include <chrono>

#ifndef IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE
#  define IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE double
#endif

namespace iporoskun::longman {

using floating_t = IPOROSKUN_LONGMAN_UNDERLYING_FLOATING_TYPE;
using duration_t = std::chrono::seconds;

struct longman_parameter {
  struct pos { // position
	double latitude; /* deg */
	double longitude; /* deg */
	double height; /* msl_orthometric_height, meter or cm*/
  } position;
  duration_t utc_offset;
};

} // namespace iporoskun::longman