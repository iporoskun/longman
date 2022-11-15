#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

#include "longman_parameter.hpp"

namespace iporoskun::longman {

struct result_type {
  double delta_g_mgal;
  double delta_g_ms2;
};

namespace details {

template<class T, std::enable_if_t<std::is_arithmetic_v<T>>...>
inline constexpr auto abs(T const& x) noexcept {
  return x < 0 ? -x : x;
}

inline constexpr double
  dms2rad(double deg, double min = 0., double sec = 0.) noexcept {
  const double angles_deg =
	details::abs(deg) + details::abs(min) / 60. + details::abs(sec) / 3600.;
  if (deg < 0. || min < 0. || sec < 0.) {
	return -angles_deg * M_PI / 180.;
  } else {
	return angles_deg * M_PI / 180.;
  }
}

inline constexpr double deg2rad(double deg) noexcept {
  return deg * M_PI / 180.;
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

inline double centuries_from_ref_date(
  std::chrono::system_clock::time_point tp,
  std::chrono::system_clock::time_point ref_date = {
	std::chrono::sys_days{ std::chrono::year{ 1899 } / 12
						   / std::chrono::day{ 31 } }
	+ std::chrono::hours{ 12 } }) noexcept {

  using std::chrono::seconds;
  using namespace std::chrono;

  const auto dp = floor<days>(tp);
  const auto diff_days = sys_days{ dp } - ref_date;

  const auto time = from_midnight<seconds>(tp);

  const auto diff = floor<seconds>(diff_days) + time;

  constexpr static auto amount_of_seconds_in_julian_century =
	100. * 365.25 * 24. * 60. * 60.;

  return static_cast<double>(diff.count())
		 / static_cast<double>(amount_of_seconds_in_julian_century);
}

} // namespace details

inline constexpr double
  dms2deg(double deg, double min = 0., double sec = 0.) noexcept {
  return details::dms2rad(deg, min, sec) * 180. / M_PI;
}

// Advances in Geophysical Methods Applied to Forensic Investigations
// https://books.google.de/books?id=NU7iDwAAQBAJ&pg=PA141&lpg=PA141&dq=Loveh2%3D0.612;+Lovek2%3D0.303;+beta%3D1%2BLoveh2-3/2*Lovek2&source=bl&ots=ZZV4OSN9e1&sig=ACfU3U1ffqTuzBzJKl4HM4VvwXFmc5Y8kA&hl=de&sa=X&ved=2ahUKEwj818LDvLftAhUi2uAKHSFOBlAQ6AEwAHoECAEQAg#v=onepage&q=Loveh2%3D0.612%3B%20Lovek2%3D0.303%3B%20beta%3D1%2BLoveh2-3%2F2*Lovek2&f=false

class longman {
  constexpr static auto a = 6.378270e8;
  constexpr static auto e_crt_2 = 0.006738;

  constexpr static auto G = 6.674e-8;
  constexpr static auto M_m = 7.3537e25;
  constexpr static auto M_s = 1.993e33;

  constexpr static auto m_s2m = 0.074804;
  constexpr static auto e_m = 0.05490;
  constexpr static auto c_m = 3.84402e10;
  constexpr static auto c_s = 1.495e13;
  constexpr static auto omega = details::dms2rad(23., 26., 21.48);
  constexpr static auto i = details::deg2rad(5.145);

  constexpr static auto love_h2 = 0.612;
  constexpr static auto love_k2 = 0.303;
  constexpr static auto beta = 1. + love_h2 - 3. / 2. * love_k2;

  constexpr static auto rev = 360.;
  constexpr static auto rev_sec = rev * 3600.;

  using time_point = std::chrono::sys_time<duration_t>;

  time_point local_time_;
  time_point utc_time_;
  duration_t utc_offset_;

  double sm_rad{};
  double pm_rad{};
  double hs_rad{};
  double N_rad{};
  double ps_rad{};
  double es{};

  longman_parameter::pos pos_deg_m;
  longman_parameter::pos pos_rad_cm;

  double r{}; // [cm] distance from observation point to the center of the Earth
  double d{}; // distance from Moon to the Earth's geocenter
  double D{}; // distance from Sun to the Earth's geocenter
  double Im_rad{}; // Inclination of the Moon's orbit to the equator
  double nu{}; // Longitude in the celestial equator of its intersection A
			   // with the Moon's orbit
  double chi_m_rad{}; // right ascension of meridian of place of observations
					  // reckoned from A
  double chi_s_rad{}; // right ascension of meridian of place of observations
					  // reckoned from the vernal equinox

public:
  longman() = default;
  explicit longman(
	const longman_parameter::pos& pos_of_msrmnt, duration_t utc_offset) {
	set_parameter(pos_of_msrmnt, utc_offset);
  }

  explicit longman(const longman_parameter& parameter) {
	set_parameter(parameter.position, parameter.utc_offset);
  }

  void set_parameter(
	const longman_parameter::pos& pos_of_msrmnt,
	duration_t utc_offset) noexcept {
	utc_offset_ = utc_offset;
	pos_deg_m = pos_of_msrmnt;
	pos_rad_cm = { details::deg2rad(pos_of_msrmnt.latitude),
				   details::deg2rad(pos_of_msrmnt.longitude),
				   pos_of_msrmnt.height * 100. };
  }

  template<class TimePoint>
  [[nodiscard]] auto operator()(const TimePoint& local_time_of_msrmnt) noexcept
	-> floating_t /*meters_per_second_squared_t*/;

  static double mean_longitude_moon(double T) noexcept;
  static double mean_longitude_lunar_perigee(double T) noexcept;
  static double mean_longitude_sun(double T) noexcept;
  static double longitude_of_moons_ascending_node(double T) noexcept;
  static double mean_longitude_solar_perigee(double T) noexcept;
  static double eccentricity_of_earths_orbit(double T) noexcept;

  static double
	distance_parameter(const longman_parameter::pos& pos_rad_cm) noexcept;
  double distance_center_moon_earth() noexcept;
  double distance_center_sun_earth() noexcept;
  double inclination_of_moon() const;
  double longitude_celestial_equator() const;

  auto get_pos_rad_cm() const noexcept { return pos_rad_cm; }

private:
  template<class TimePoint>
  void set_time(const TimePoint& local_time) noexcept;

  void calc_longitude_and_eccentricity(double T) noexcept;

  auto calculate_acceleration() -> floating_t
	/*-> meters_per_second_squared_t*/;
};

template<class TimePoint>
inline auto longman::operator()(const TimePoint& local_time_of_msrmnt) noexcept
  -> floating_t /*meters_per_second_squared_t*/
{
  using namespace std::chrono;

  set_time<TimePoint>(local_time_of_msrmnt);
  const auto T =
	details::centuries_from_ref_date(floor<seconds>(local_time_of_msrmnt));

  calc_longitude_and_eccentricity(T);
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