#pragma once
#include <cmath>
#include <numbers>
#include <concepts>

namespace iporoskun::longman {

namespace detail {
// Implementation from
// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/

template<class T, class /*TagType*/>
class named_type {
public:
  explicit named_type(T const& value) noexcept
	: value_(value) {}
  explicit named_type(T&& value) noexcept
	: value_(std::move(value)) {}

  named_type() noexcept = default;
  named_type(const named_type&) noexcept = default;
  named_type(named_type&&) noexcept = default;
  named_type& operator=(named_type const&) noexcept = default;
  named_type& operator=(named_type&&) noexcept = default;
  operator T() const noexcept { return value_; };

  T& get() noexcept { return value_; }
  T const& get() const noexcept { return value_; }

private:
  T value_;
};

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

namespace detail {

template<std::floating_point T>
inline constexpr auto abs(T const& x) noexcept {
  return x < static_cast<T>(0) ? -x : x;
}

template<std::floating_point floating_t>
inline constexpr floating_t dms_to_rad(
  floating_t deg, floating_t min = 0., floating_t sec = 0.) noexcept {
  const floating_t angles_deg =
	detail::abs(deg) + detail::abs(min) / 60. + detail::abs(sec) / 3600.;
  if (deg < 0. || min < 0. || sec < 0.) {
	return -angles_deg
		   * std::numbers::pi_v<floating_t> / static_cast<floating_t>(180);
  } else {
	return angles_deg
		   * std::numbers::pi_v<floating_t> / static_cast<floating_t>(180);
  }
}

template<std::floating_point floating_t>
inline constexpr auto deg_to_rad(floating_t deg) noexcept -> floating_t {
  return deg * std::numbers::pi_v<floating_t> / static_cast<floating_t>(180);
}

template<std::floating_point floating_t>
inline constexpr auto rad_to_deg(floating_t rad) noexcept -> floating_t {
  return rad / std::numbers::pi_v<floating_t> * static_cast<floating_t>(180);
}

template<class TimePointType>
inline auto time_since_midnight(TimePointType tp) noexcept {
  const auto dp = std::chrono::floor<std::chrono::days>(tp);
  return std::chrono::hh_mm_ss(tp - dp);
}

template<class Duration, class TimePointType>
inline auto from_midnight(TimePointType tp) noexcept {
  const auto time = time_since_midnight(std::move(tp));
  return std::chrono::floor<Duration>(time.to_duration());
}

template<std::floating_point floating_t, class TimePointType>
inline floating_t
  julian_centuries_from_reference_date(TimePointType tp) noexcept {

  using namespace std::chrono;

  static constexpr system_clock::time_point reference_date = {
	sys_days{ year{ 1899 } / 12 / day{ 31 } } + hours{ 12 }
  };

  const auto dp = floor<days>(tp);
  const auto diff_days = sys_days{ dp } - reference_date;
  const auto time_diff = from_midnight<std::chrono::seconds>(tp);
  const auto total_diff = floor<std::chrono::seconds>(diff_days) + time_diff;

  constexpr static auto number_of_seconds_in_julian_century =
	static_cast<floating_t>(3155760000); /*100*365.25*24*60*60*/

  return static_cast<floating_t>(total_diff.count())
		 / number_of_seconds_in_julian_century;
}

template<std::floating_point floating_t>
inline constexpr floating_t degree_minute_second_to_degree(
  floating_t deg, floating_t min, floating_t sec) noexcept {
  return detail::dms_to_rad(deg, min, sec) * 180.
		 / std::numbers::pi_v<floating_t>;
}


auto mean_longitude_moon(std::floating_point auto time) noexcept;
auto mean_longitude_lunar_perigee(std::floating_point auto time) noexcept;
auto mean_longitude_sun(std::floating_point auto time) noexcept;
auto longitude_of_moons_ascending_node(std::floating_point auto time) noexcept;
auto mean_longitude_solar_perigee(std::floating_point auto time) noexcept;
auto eccentricity_of_earths_orbit(std::floating_point auto time) noexcept;

template<std::floating_point FloatingType>
inline auto position_from_deg_meter_to_rad_cm(position<FloatingType> const& pos)
  -> position<FloatingType> {
  return position<FloatingType>{
	latitude<FloatingType>{ detail::deg_to_rad<FloatingType>(pos.latitude) },
	longitude<FloatingType>{ detail::deg_to_rad<FloatingType>(pos.longitude) },
	height<FloatingType>{ pos.height * 100. }
  };
}

template<
  std::floating_point FloatingType = double,
  class DurationType = std::chrono::seconds>
class longman_impl {
  using time_point_t = std::chrono::system_clock::time_point;
  using position_t = position<FloatingType>;
  using floating_t = FloatingType;

public:
  explicit longman_impl(const position_t& position) {
	pos_rad_cm = position_from_deg_meter_to_rad_cm<floating_t>(position);
  }

  [[nodiscard]] auto operator()(const time_point_t& utc_time) noexcept
	-> floating_t /*meters_per_second_squared_t*/;

  static floating_t
	distance_to_earth_centre(const position_t& pos_rad_cm) noexcept;
  floating_t distance_center_moon_earth() const noexcept;
  floating_t distance_center_sun_earth() const noexcept;
  floating_t inclination_of_moon() const;
  floating_t longitude_celestial_equator() const;

private:
  void calc_longitude_and_eccentricity(floating_t time) noexcept;

  auto calculate_acceleration(const time_point_t& utc_time) const -> floating_t
	/*-> meters_per_second_squared_t*/;

  position_t pos_rad_cm;

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
};

template<std::floating_point FloatingType, class Duration>
inline auto longman_impl<FloatingType, Duration>::operator()(
  const time_point_t& utc_time) noexcept
  -> floating_t /*meters_per_second_squared_t*/
{
  using namespace std::chrono;
  const auto time = julian_centuries_from_reference_date<FloatingType>(
	floor<Duration>(utc_time));

  calc_longitude_and_eccentricity(time);
  r = distance_to_earth_centre(pos_rad_cm);
  d = distance_center_moon_earth();
  D = distance_center_sun_earth();
  Im_rad = inclination_of_moon();

  return calculate_acceleration(utc_time);
}

// Implementation: previously longman.cpp

template<std::floating_point floating_t = double>
struct constant {
  // Advances in Geophysical Methods Applied to Forensic Investigations
  // https://shorturl.at/azJ49
  static constexpr floating_t a = 6.378270e8;
  static constexpr floating_t e_crt_2 = 0.006738;

  static constexpr floating_t G = 6.674e-8;
  static constexpr floating_t M_m = 7.3537e25;
  static constexpr floating_t M_s = 1.993e33;

  static constexpr floating_t m_s2m = 0.074804;
  static constexpr floating_t e_m = 0.05490;
  static constexpr floating_t c_m = 3.84402e10;
  static constexpr floating_t c_s = 1.495e13;
  static constexpr floating_t omega = detail::dms_to_rad(23., 26., 21.48);
  static constexpr floating_t i = detail::deg_to_rad(5.145);

  static constexpr floating_t love_h2 = 0.612;
  static constexpr floating_t love_k2 = 0.303;
  static constexpr floating_t beta = 1. + love_h2 - 3. / 2. * love_k2;

  static constexpr floating_t rev_sec = 360. * 3600.;
};

template<std::floating_point FloatingType = double>
auto mean_longitude_moon(FloatingType time) noexcept { // sm_rad
  // const auto sm_rad = detail::dms_to_rad(270., 26., 14.72) +
  // detail::deg_to_rad((1336. * rev_sec + 1108411.20) / 3600.)*time +
  // detail::deg_to_rad(9.09 / 3600.)*std::pow(time, 2.) +
  // detail::deg_to_rad(0.0068 / 3600.)*std::pow(time, 3.);

  const auto sm_rad =
	detail::dms_to_rad(270., 26., 11.72)
	+ detail::deg_to_rad(
		(1336. * constant<FloatingType>::rev_sec + 1108406.05) / 3600.)
		* time
	+ detail::deg_to_rad(7.128 / 3600.) * std::pow(time, 2.)
	+ detail::deg_to_rad(0.0072 / 3600.) * std::pow(time, 3.);

  return sm_rad;
}

template<std::floating_point FloatingType = double>
auto mean_longitude_lunar_perigee(FloatingType time) noexcept { // pm_rad
  // const auto pm_rad = detail::dms_to_rad(334., 19., 40.87) +
  // detail::deg_to_rad((11. *  constants::rev_sec + 392515.94) / 3600.)*time -
  // detail::deg_to_rad(37.24 / 3600.)*std::pow(time, 2.) -
  // detail::deg_to_rad(0.045 / 3600.)*std::pow(time, 3.);

  const auto pm_rad =
	detail::dms_to_rad(334., 19., 46.42)
	+ detail::deg_to_rad(
		(11. * constant<FloatingType>::rev_sec + 392522.51) / 3600.)
		* time
	- detail::deg_to_rad(37.15 / 3600.) * std::pow(time, 2.)
	- detail::deg_to_rad(0.036 / 3600.) * std::pow(time, 3.);

  return pm_rad;
}

auto mean_longitude_sun(std::floating_point auto time) noexcept { // hs_rad
  // const auto hs_rad = detail::dms_to_rad(279., 41., 48.04) +
  // detail::deg_to_rad(129602768.13 / 3600.)*time + detail::deg_to_rad(1.089
  // / 3600.)*std::pow(time, 2.);

  const auto hs_rad = detail::dms_to_rad(279., 41., 48.05)
					  + detail::deg_to_rad(129602768.11 / 3600.) * time
					  + detail::deg_to_rad(1.080 / 3600.) * std::pow(time, 2.);

  return hs_rad;
}

template<std::floating_point FloatingType = double>
auto longitude_of_moons_ascending_node(
  std::floating_point auto time) noexcept { // N_rad
  // const auto N_rad = detail::dms_to_rad(259., 10., 57.12) -
  // detail::deg_to_rad((5. * rev_sec + 482912.63) / 3600.)*time +
  // detail::deg_to_rad(7.58 / 3600.)*std::pow(time, 2.) +
  // detail::deg_to_rad(0.008 / 3600.)*std::pow(time, 3.);

  const auto N_rad =
	detail::dms_to_rad(259., 10., 59.81)
	- detail::deg_to_rad(
		(5. * constant<FloatingType>::rev_sec + 482911.24) / 3600.)
		* time
	+ detail::deg_to_rad(7.48 / 3600.) * std::pow(time, 2.)
	+ detail::deg_to_rad(0.007 / 3600.) * std::pow(time, 3.);
  return N_rad;
}

template<std::floating_point FloatingType = double>
auto mean_longitude_solar_perigee(FloatingType time) noexcept { // ps_rad
  // const auto ps_rad = detail::dms_to_rad(281., 13., 15.00) +
  // detail::deg_to_rad(6189.03 / 3600.)*time + detail::deg_to_rad(1.63 /
  // 3600.)*std::pow(time, 2.) + detail::deg_to_rad(0.012 / 3600.)*
  // std::pow(time, 3.);

  const auto ps_rad = detail::dms_to_rad(281., 13., 14.99)
					  + detail::deg_to_rad(6188.47 / 3600.) * time
					  + detail::deg_to_rad(1.62 / 3600.) * std::pow(time, 2.)
					  + detail::deg_to_rad(0.011 / 3600.) * std::pow(time, 3.);

  return ps_rad;
}

auto eccentricity_of_earths_orbit(
  std::floating_point auto time) noexcept { // es
  return 0.01675104 - 0.00004180 * time - 0.000000126 * std::pow(time, 2.);
}

template<std::floating_point F, class D>
void longman_impl<F, D>::calc_longitude_and_eccentricity(
  floating_t time) noexcept {
  sm_rad = detail::mean_longitude_moon<F>(time);
  pm_rad = detail::mean_longitude_lunar_perigee<F>(time);
  hs_rad = detail::mean_longitude_sun<F>(time);
  N_rad = detail::longitude_of_moons_ascending_node<F>(time);
  ps_rad = detail::mean_longitude_solar_perigee<F>(time);
  es = detail::eccentricity_of_earths_orbit<F>(time);
}

template<std::floating_point FloatingType, class D>
auto longman_impl<FloatingType, D>::calculate_acceleration(
  const time_point_t& utc_time) const
  -> FloatingType /* meters_per_second_squared_t*/
{
  using constant = constant<FloatingType>;
  using fhours_t = std::chrono::duration<floating_t, std::ratio<3600>>;
  const auto t0 = detail::from_midnight<fhours_t>(utc_time).count();
  const auto t_rad = detail::deg_to_rad<FloatingType>(
	15. * (t0 - 12.) + detail::rad_to_deg<FloatingType>(pos_rad_cm.longitude));

  // Longitude in the celestial equator of its intersection A with the Moon's
  // orbit
  const auto nu = longitude_celestial_equator();

  const auto cos_alpha =
	cos(N_rad) * cos(nu) + sin(N_rad) * sin(nu) * cos(constant::omega);
  const auto sin_alpha = (sin(constant::omega) * sin(N_rad)) / sin(Im_rad);
  const auto alpha = 2 * atan(sin_alpha / (1. + cos_alpha));

  const auto xi = N_rad - alpha;
  const auto sigma = sm_rad - xi;

  const auto L_m_rad =
	sigma + 2. * constant::e_m * sin(sm_rad - pm_rad)
	+ 5. / 4. * pow(constant::e_m, 2.) * sin(2. * (sm_rad - pm_rad))
	+ 15. / 4. * constant::m_s2m * constant::e_m
		* sin(sm_rad - 2. * hs_rad + pm_rad)
	+ 11. / 8. * pow(constant::m_s2m, 2.) * sin(2. * (sm_rad - hs_rad));

  const auto L_sun_rad = hs_rad + 2. * es * sin(hs_rad - ps_rad);

  // Zenith angle of the Moon
  // e.g. right ascension of meridian of place of observations reckoned from A
  const auto chi_m_rad = t_rad + hs_rad - nu;
  const auto cos_Zm =
	sin(pos_rad_cm.latitude) * sin(Im_rad) * sin(L_m_rad)
	+ cos(pos_rad_cm.latitude)
		* (pow(cos(Im_rad / 2.), 2.) * cos(L_m_rad - chi_m_rad) + pow(sin(Im_rad / 2.), 2.) * cos(L_m_rad + chi_m_rad));

  // Zenith angle of the Sun
  // e.g. right ascension of meridian of place of observations reckoned from the
  // vernal equinox
  const auto chi_s_rad = t_rad + hs_rad;
  const auto cos_Zs =
	sin(pos_rad_cm.latitude) * sin(constant::omega) * sin(L_sun_rad)
	+ cos(pos_rad_cm.latitude)
		* (pow(cos(constant::omega / 2.), 2.) * cos(L_sun_rad - chi_s_rad) + pow(sin(constant::omega / 2.), 2.) * cos(L_sun_rad + chi_s_rad));

  // vertical component of tidal acceleration due to the Moon
  const auto gm = constant::G * constant::M_m * r * pow(1. / d, 3.)
					* (3. * pow(cos_Zm, 2.) - 1.)
				  + 1.5 * constant::G * constant::M_m * pow(r, 2.)
					  * pow(1. / d, 4.) * (5. * pow(cos_Zm, 3.) - 3. * cos_Zm);

  // vertical component of tidal acceleration due to the Sun :
  const auto gs = constant::G * constant::M_s * r * pow(1. / D, 3.)
				  * (3. * pow(cos_Zs, 2.) - 1.);

  const auto g0_gal = (gm + gs) * constant::beta;
  // return { g0_gal * 1000., g0_gal / 100. };
  return /*meters_per_second_squared_t*/ floating_t{ g0_gal / 100. };
}

template<std::floating_point F, class D>
auto longman_impl<F, D>::distance_to_earth_centre(
  const position_t& pos_rad_cm) noexcept -> floating_t { // r
  const auto C_2 =
	1. / (1. + constant<F>::e_crt_2 * pow(sin(pos_rad_cm.latitude.get()), 2.));
  const auto C = sqrt(C_2);
  return C * constant<F>::a + pos_rad_cm.height.get();
}

template<std::floating_point F, class D>
auto longman_impl<F, D>::distance_center_moon_earth() const noexcept
  -> floating_t { // d
  using constant = constant<F>;
  const auto a_moon = 1. / (constant::c_m * (1. - pow(constant::e_m, 2.)));
  const auto recip_d =
	1. / constant::c_m + a_moon * constant::e_m * cos(sm_rad - pm_rad)
	+ a_moon * pow(constant::e_m, 2.) * cos(2. * (sm_rad - pm_rad))
	+ (15. / 8.) * a_moon * constant::m_s2m * constant::e_m
		* cos(sm_rad - 2. * hs_rad + pm_rad)
	+ a_moon * pow(constant::m_s2m, 2.) * cos(2. * (sm_rad - hs_rad));
  return 1. / recip_d;
}

template<std::floating_point F, class D>
auto longman_impl<F, D>::distance_center_sun_earth() const noexcept
  -> floating_t { // D
  using constant = constant<F>;
  const auto a_sun = 1. / (constant::c_s * (1. - pow(es, 2.)));
  const auto recip_D =
	((1. / constant::c_s) + a_sun * es * cos(hs_rad - ps_rad));
  return static_cast<floating_t>(1) / recip_D;
}

template<std::floating_point F, class D>
auto longman_impl<F, D>::inclination_of_moon() const -> floating_t { // Im_rad
  using constant = constant<F>;
  const auto inc_moon_rad = acos(
	cos(constant::omega) * cos(constant::i)
	- sin(constant::omega) * sin(constant::i) * cos(N_rad));
#if not defined(IPOROSKUN_LONGMAN_DISABLE_RUNTIME_CHECKS)
  if (
	((inc_moon_rad * 180 / std::numbers::pi) < 18)
	|| ((inc_moon_rad * 180 / std::numbers::pi) > 28)) {
	throw std::runtime_error("Inclination of the Moon's orbit to the equator "
							 "is outside of the normal range.");
  }
#endif
  return inc_moon_rad;
}

template<std::floating_point F, class D>
auto longman_impl<F, D>::longitude_celestial_equator() const
  -> floating_t { // nu
  const auto nu_t = asin(sin(constant<F>::i) * sin(N_rad) / sin(Im_rad));
#if not defined(IPOROSKUN_LONGMAN_DISABLE_RUNTIME_CHECKS)
  if (
	((nu_t * 180 / std::numbers::pi_v<floating_t>) < -15)
	|| ((nu_t * 180 / std::numbers::pi_v<floating_t>) > 15)) {
	throw std::runtime_error(
	  "Longitude in the celestial equator is not correct.");
  }
#endif
  return nu_t;
}

// End implementation

} // namespace detail

template<
  std::floating_point FloatingType = double,
  class DurationType = std::chrono::seconds,
  class TimePoint = std::chrono::system_clock::time_point>
[[nodiscard]] inline auto
  longman(const position<FloatingType>& position, const TimePoint& utc_time)
	-> FloatingType {
  return detail::longman_impl<FloatingType, DurationType>(position)(utc_time);
}

} // namespace iporoskun::longman