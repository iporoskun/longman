#pragma once
#include <cmath>
#include <numbers>

#include "longman_parameter.hpp"

namespace iporoskun::longman {

struct result_type {
  floating_t delta_g_mgal;
  floating_t delta_g_ms2;
};

namespace details {

using namespace std::numbers;

template<std::floating_point T>
inline constexpr auto abs(T const& x) noexcept {
  return x < 0 ? -x : x;
}

inline constexpr floating_t dms_to_rad(
  floating_t deg, floating_t min = 0., floating_t sec = 0.) noexcept {
  const floating_t angles_deg =
	details::abs(deg) + details::abs(min) / 60. + details::abs(sec) / 3600.;
  if (deg < 0. || min < 0. || sec < 0.) {
	return -angles_deg * pi_v<floating_t> / static_cast<floating_t>(180);
  } else {
	return angles_deg * pi_v<floating_t> / static_cast<floating_t>(180);
  }
}

inline constexpr floating_t deg_to_rad(floating_t deg) noexcept {
  return deg * pi_v<floating_t> / static_cast<floating_t>(180);
}

inline constexpr floating_t rad2deg(floating_t rad) noexcept {
  return rad / pi_v<floating_t> * static_cast<floating_t>(180);
}

inline auto
  time_since_midnight(std::chrono::system_clock::time_point tp) noexcept {
  const auto dp = std::chrono::floor<std::chrono::days>(tp);
  return std::chrono::hh_mm_ss(tp - dp);
}

template<class Duration>
inline auto from_midnight(std::chrono::system_clock::time_point tp) noexcept {
  const auto time = time_since_midnight(std::move(tp));
  return std::chrono::floor<Duration>(time.to_duration());
}

inline floating_t julian_centuries_from_reference_date(
  std::chrono::system_clock::time_point tp) noexcept {

  using std::chrono::seconds;
  using namespace std::chrono;

  static constexpr system_clock::time_point reference_date = {
	sys_days{ year{ 1899 } / 12 / day{ 31 } } + hours{ 12 }
  };

  const auto dp = floor<days>(tp);
  const auto diff_days = sys_days{ dp } - reference_date;
  const auto time_diff = from_midnight<seconds>(tp);
  const auto total_diff = floor<seconds>(diff_days) + time_diff;

  constexpr static auto number_of_seconds_in_julian_century =
	static_cast<floating_t>(3155760000); /*100*365.25*24*60*60*/

  return static_cast<floating_t>(total_diff.count())
		 / number_of_seconds_in_julian_century;
}

} // namespace details

inline constexpr floating_t degree_minute_second_to_degree(
  floating_t deg, floating_t min, floating_t sec) noexcept {
  return details::dms_to_rad(deg, min, sec) * 180.
		 / std::numbers::pi_v<floating_t>;
}

namespace constants {
// Advances in Geophysical Methods Applied to Forensic Investigations
// https://shorturl.at/azJ49
inline constexpr floating_t a = 6.378270e8;
inline constexpr floating_t e_crt_2 = 0.006738;

inline constexpr floating_t G = 6.674e-8;
inline constexpr floating_t M_m = 7.3537e25;
inline constexpr floating_t M_s = 1.993e33;

inline constexpr floating_t m_s2m = 0.074804;
inline constexpr floating_t e_m = 0.05490;
inline constexpr floating_t c_m = 3.84402e10;
inline constexpr floating_t c_s = 1.495e13;
inline constexpr floating_t omega = details::dms_to_rad(23., 26., 21.48);
inline constexpr floating_t i = details::deg_to_rad(5.145);

inline constexpr floating_t love_h2 = 0.612;
inline constexpr floating_t love_k2 = 0.303;
inline constexpr floating_t beta = 1. + love_h2 - 3. / 2. * love_k2;

inline constexpr floating_t rev_sec = 360. * 3600.;
} // namespace constants

namespace details {
floating_t mean_longitude_moon(floating_t time) noexcept;
floating_t mean_longitude_lunar_perigee(floating_t time) noexcept;
floating_t mean_longitude_sun(floating_t time) noexcept;
floating_t longitude_of_moons_ascending_node(floating_t time) noexcept;
floating_t mean_longitude_solar_perigee(floating_t time) noexcept;
floating_t eccentricity_of_earths_orbit(floating_t time) noexcept;

inline longman_parameter::position_t
  position_from_deg_meter_to_rad_cm(longman_parameter::position_t const& pos) {
  return longman_parameter::position_t{
	latitude{ details::deg_to_rad(pos.latitude) },
	longitude{ details::deg_to_rad(pos.longitude) }, height{ pos.height * 100. }
  };
}
} // namespace details


class longman {
  using time_point = std::chrono::system_clock::time_point;

public:
  longman(longman_parameter::position_t const& measurement_position) {
	pos_rad_cm =
	  details::position_from_deg_meter_to_rad_cm(measurement_position);
  }

  explicit longman(const longman_parameter& parameter)
	: longman(parameter.position) {}

  [[nodiscard]] auto operator()(const time_point& utc_time) noexcept
	-> floating_t /*meters_per_second_squared_t*/;

  static floating_t distance_to_earth_centre(
	const longman_parameter::position_t& pos_rad_cm) noexcept;
  floating_t distance_center_moon_earth() noexcept;
  floating_t distance_center_sun_earth() noexcept;
  floating_t inclination_of_moon() const;
  floating_t longitude_celestial_equator() const;

private:
  void calc_longitude_and_eccentricity(floating_t time) noexcept;

  auto calculate_acceleration(const time_point& utc_time) -> floating_t
	/*-> meters_per_second_squared_t*/;

  longman_parameter::position_t pos_rad_cm;

  floating_t sm_rad;
  floating_t pm_rad;
  floating_t hs_rad;
  floating_t N_rad;
  floating_t ps_rad;
  floating_t es;

  floating_t
	r; // [cm] distance from observation point to the center of the Earth
  floating_t d; // distance from Moon to the Earth's geocenter
  floating_t D; // distance from Sun to the Earth's geocenter
  floating_t Im_rad; // Inclination of the Moon's orbit to the equator
  floating_t nu; // Longitude in the celestial equator of its
				 // intersection A with the Moon's orbit
  floating_t chi_m_rad; // right ascension of meridian of place of
						// observations reckoned from A
  floating_t chi_s_rad; // right ascension of meridian of place of
						// observations reckoned from the vernal equinox
};

inline auto longman::operator()(const time_point& utc_time) noexcept
  -> floating_t /*meters_per_second_squared_t*/
{
  using namespace std::chrono;
  const auto time =
	details::julian_centuries_from_reference_date(floor<seconds>(utc_time));

  calc_longitude_and_eccentricity(time);
  r = distance_to_earth_centre(pos_rad_cm);
  d = distance_center_moon_earth();
  D = distance_center_sun_earth();
  Im_rad = inclination_of_moon();
  nu = longitude_celestial_equator();

  return calculate_acceleration(utc_time);
}

} // namespace iporoskun::longman