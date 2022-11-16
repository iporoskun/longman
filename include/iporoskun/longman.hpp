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

inline constexpr floating_t
  dms2rad(floating_t deg, floating_t min = 0., floating_t sec = 0.) noexcept {
  const floating_t angles_deg =
	details::abs(deg) + details::abs(min) / 60. + details::abs(sec) / 3600.;
  if (deg < 0. || min < 0. || sec < 0.) {
	return -angles_deg * pi_v<floating_t> / static_cast<floating_t>(180);
  } else {
	return angles_deg * pi_v<floating_t> / static_cast<floating_t>(180);
  }
}

inline constexpr floating_t deg2rad(floating_t deg) noexcept {
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

inline floating_t julian_centuries_from_ref_date(
  std::chrono::system_clock::time_point tp,
  std::chrono::system_clock::time_point ref_date = {
	std::chrono::sys_days{ std::chrono::year{ 1899 } / 12
						   / std::chrono::day{ 31 } }
	+ std::chrono::hours{ 12 } }) noexcept {

  using std::chrono::seconds;
  using namespace std::chrono;

  const auto dp = floor<days>(tp);
  const auto diff_days = sys_days{ dp } - ref_date;
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
  return details::dms2rad(deg, min, sec) * 180.
		 / std::numbers::pi_v<floating_t>;
}

// Advances in Geophysical Methods Applied to Forensic Investigations
// https://books.google.de/books?id=NU7iDwAAQBAJ&pg=PA141&lpg=PA141&dq=Loveh2%3D0.612;+Lovek2%3D0.303;+beta%3D1%2BLoveh2-3/2*Lovek2&source=bl&ots=ZZV4OSN9e1&sig=ACfU3U1ffqTuzBzJKl4HM4VvwXFmc5Y8kA&hl=de&sa=X&ved=2ahUKEwj818LDvLftAhUi2uAKHSFOBlAQ6AEwAHoECAEQAg#v=onepage&q=Loveh2%3D0.612%3B%20Lovek2%3D0.303%3B%20beta%3D1%2BLoveh2-3%2F2*Lovek2&f=false


class longman {
public:
  longman() = default;
  explicit longman(
	const longman_parameter::position_t& pos_of_msrmnt, duration_t utc_offset) {
	set_parameter(pos_of_msrmnt, utc_offset);
  }

  explicit longman(const longman_parameter& parameter) {
	set_parameter(parameter.position, parameter.utc_offset);
  }

  void set_parameter(
	const longman_parameter::position_t& pos_of_msrmnt,
	duration_t utc_offset) noexcept {
	utc_offset_ = utc_offset;
	pos_rad_cm = longman_parameter::position_t{
	  latitude{ details::deg2rad(pos_of_msrmnt.latitude) },
	  longitude{ details::deg2rad(pos_of_msrmnt.longitude) },
	  height{ pos_of_msrmnt.height * 100. }
	};
  }

  template<class TimePoint>
  [[nodiscard]] auto operator()(const TimePoint& local_time_of_msrmnt) noexcept
	-> floating_t /*meters_per_second_squared_t*/;

  static floating_t mean_longitude_moon(floating_t T) noexcept;
  static floating_t mean_longitude_lunar_perigee(floating_t T) noexcept;
  static floating_t mean_longitude_sun(floating_t T) noexcept;
  static floating_t longitude_of_moons_ascending_node(floating_t T) noexcept;
  static floating_t mean_longitude_solar_perigee(floating_t T) noexcept;
  static floating_t eccentricity_of_earths_orbit(floating_t T) noexcept;

  static floating_t distance_parameter(
	const longman_parameter::position_t& pos_rad_cm) noexcept;
  floating_t distance_center_moon_earth() noexcept;
  floating_t distance_center_sun_earth() noexcept;
  floating_t inclination_of_moon() const;
  floating_t longitude_celestial_equator() const;

  auto get_pos_rad_cm() const noexcept { return pos_rad_cm; }

private:
  template<class TimePoint>
  void set_time(const TimePoint& local_time) noexcept;

  void calc_longitude_and_eccentricity(floating_t T) noexcept;

  auto calculate_acceleration() -> floating_t
	/*-> meters_per_second_squared_t*/;

  constexpr static floating_t a = 6.378270e8;
  constexpr static floating_t e_crt_2 = 0.006738;

  constexpr static floating_t G = 6.674e-8;
  constexpr static floating_t M_m = 7.3537e25;
  constexpr static floating_t M_s = 1.993e33;

  constexpr static floating_t m_s2m = 0.074804;
  constexpr static floating_t e_m = 0.05490;
  constexpr static floating_t c_m = 3.84402e10;
  constexpr static floating_t c_s = 1.495e13;
  constexpr static floating_t omega = details::dms2rad(23., 26., 21.48);
  constexpr static floating_t i = details::deg2rad(5.145);

  constexpr static floating_t love_h2 = 0.612;
  constexpr static floating_t love_k2 = 0.303;
  constexpr static floating_t beta = 1. + love_h2 - 3. / 2. * love_k2;

  constexpr static floating_t rev_sec = 360. * 3600.;

  using time_point = std::chrono::sys_time<duration_t>;
  time_point local_time_;
  time_point utc_time_;
  duration_t utc_offset_;

  floating_t sm_rad;
  floating_t pm_rad;
  floating_t hs_rad;
  floating_t N_rad;
  floating_t ps_rad;
  floating_t es;

  longman_parameter::position_t pos_rad_cm;

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

template<class TimePoint>
inline auto longman::operator()(const TimePoint& local_time_of_msrmnt) noexcept
  -> floating_t /*meters_per_second_squared_t*/
{
  using namespace std::chrono;

  set_time<TimePoint>(local_time_of_msrmnt);
  const auto time = details::julian_centuries_from_ref_date(
	floor<seconds>(local_time_of_msrmnt));

  calc_longitude_and_eccentricity(time);
  r = distance_parameter(pos_rad_cm);
  d = distance_center_moon_earth();
  D = distance_center_sun_earth();
  Im_rad = inclination_of_moon();
  nu = longitude_celestial_equator();

  return calculate_acceleration();
}

template<class TimePoint>
inline void longman::set_time(const TimePoint& local_time) noexcept {
  local_time_ = std::chrono::floor<std::chrono::seconds>(local_time);
  utc_time_ = local_time_ - utc_offset_;
}

} // namespace iporoskun::longman