#include <exception>
#include "iporoskun/longman.hpp"

namespace iporoskun::longman {

namespace constant {
// Advances in Geophysical Methods Applied to Forensic Investigations
// https://shorturl.at/azJ49
constexpr floating_t a = 6.378270e8;
constexpr floating_t e_crt_2 = 0.006738;

constexpr floating_t G = 6.674e-8;
constexpr floating_t M_m = 7.3537e25;
constexpr floating_t M_s = 1.993e33;

constexpr floating_t m_s2m = 0.074804;
constexpr floating_t e_m = 0.05490;
constexpr floating_t c_m = 3.84402e10;
constexpr floating_t c_s = 1.495e13;
constexpr floating_t omega = details::dms_to_rad(23., 26., 21.48);
constexpr floating_t i = details::deg_to_rad(5.145);

constexpr floating_t love_h2 = 0.612;
constexpr floating_t love_k2 = 0.303;
constexpr floating_t beta = 1. + love_h2 - 3. / 2. * love_k2;

constexpr floating_t rev_sec = 360. * 3600.;
} // namespace constant

floating_t details::mean_longitude_moon(floating_t time) noexcept { // sm_rad
  // const auto sm_rad = details::dms_to_rad(270., 26., 14.72) +
  // details::deg_to_rad((1336. * rev_sec + 1108411.20) / 3600.)*time +
  // details::deg_to_rad(9.09 / 3600.)*std::pow(time, 2.) +
  // details::deg_to_rad(0.0068 / 3600.)*std::pow(time, 3.);

  const auto sm_rad =
	details::dms_to_rad(270., 26., 11.72)
	+ details::deg_to_rad((1336. * constant::rev_sec + 1108406.05) / 3600.)
		* time
	+ details::deg_to_rad(7.128 / 3600.) * std::pow(time, 2.)
	+ details::deg_to_rad(0.0072 / 3600.) * std::pow(time, 3.);

  return sm_rad;
}

floating_t
  details::mean_longitude_lunar_perigee(floating_t time) noexcept { // pm_rad
  // const auto pm_rad = details::dms_to_rad(334., 19., 40.87) +
  // details::deg_to_rad((11. *  constants::rev_sec + 392515.94) / 3600.)*time -
  // details::deg_to_rad(37.24 / 3600.)*std::pow(time, 2.) -
  // details::deg_to_rad(0.045 / 3600.)*std::pow(time, 3.);

  const auto pm_rad =
	details::dms_to_rad(334., 19., 46.42)
	+ details::deg_to_rad((11. * constant::rev_sec + 392522.51) / 3600.) * time
	- details::deg_to_rad(37.15 / 3600.) * std::pow(time, 2.)
	- details::deg_to_rad(0.036 / 3600.) * std::pow(time, 3.);

  return pm_rad;
}

floating_t details::mean_longitude_sun(floating_t time) noexcept { // hs_rad
  // const auto hs_rad = details::dms_to_rad(279., 41., 48.04) +
  // details::deg_to_rad(129602768.13 / 3600.)*time + details::deg_to_rad(1.089
  // / 3600.)*std::pow(time, 2.);

  const auto hs_rad = details::dms_to_rad(279., 41., 48.05)
					  + details::deg_to_rad(129602768.11 / 3600.) * time
					  + details::deg_to_rad(1.080 / 3600.) * std::pow(time, 2.);

  return hs_rad;
}

floating_t details::longitude_of_moons_ascending_node(
  floating_t time) noexcept { // N_rad
  // const auto N_rad = details::dms_to_rad(259., 10., 57.12) -
  // details::deg_to_rad((5. * rev_sec + 482912.63) / 3600.)*time +
  // details::deg_to_rad(7.58 / 3600.)*std::pow(time, 2.) +
  // details::deg_to_rad(0.008 / 3600.)*std::pow(time, 3.);

  const auto N_rad =
	details::dms_to_rad(259., 10., 59.81)
	- details::deg_to_rad((5. * constant::rev_sec + 482911.24) / 3600.) * time
	+ details::deg_to_rad(7.48 / 3600.) * std::pow(time, 2.)
	+ details::deg_to_rad(0.007 / 3600.) * std::pow(time, 3.);
  return N_rad;
}

floating_t
  details::mean_longitude_solar_perigee(floating_t time) noexcept { // ps_rad
  // const auto ps_rad = details::dms_to_rad(281., 13., 15.00) +
  // details::deg_to_rad(6189.03 / 3600.)*time + details::deg_to_rad(1.63 /
  // 3600.)*std::pow(time, 2.) + details::deg_to_rad(0.012 / 3600.)*
  // std::pow(time, 3.);

  const auto ps_rad = details::dms_to_rad(281., 13., 14.99)
					  + details::deg_to_rad(6188.47 / 3600.) * time
					  + details::deg_to_rad(1.62 / 3600.) * std::pow(time, 2.)
					  + details::deg_to_rad(0.011 / 3600.) * std::pow(time, 3.);

  return ps_rad;
}

floating_t
  details::eccentricity_of_earths_orbit(floating_t time) noexcept { // es
  return 0.01675104 - 0.00004180 * time - 0.000000126 * std::pow(time, 2.);
}

void longman::calc_longitude_and_eccentricity(floating_t time) noexcept {
  sm_rad = details::mean_longitude_moon(time);
  pm_rad = details::mean_longitude_lunar_perigee(time);
  hs_rad = details::mean_longitude_sun(time);
  N_rad = details::longitude_of_moons_ascending_node(time);
  ps_rad = details::mean_longitude_solar_perigee(time);
  es = details::eccentricity_of_earths_orbit(time);
}

auto longman::calculate_acceleration(const time_point& utc_time) const
  -> floating_t /* meters_per_second_squared_t*/
{
  using fhours_t = std::chrono::duration<floating_t, std::ratio<3600>>;
  const auto t0 = details::from_midnight<fhours_t>(utc_time).count();
  const auto t_rad = details::deg_to_rad(
	15. * (t0 - 12.) + details::rad_to_deg(pos_rad_cm.longitude));

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


floating_t
  longman::distance_to_earth_centre(const position& pos_rad_cm) noexcept { // r
  const auto C_2 =
	1. / (1. + constant::e_crt_2 * pow(sin(pos_rad_cm.latitude.get()), 2.));
  const auto C = sqrt(C_2);
  return C * constant::a + pos_rad_cm.height.get();
}

floating_t longman::distance_center_moon_earth() const noexcept { // d
  using namespace constant;
  const auto a_moon = 1. / (c_m * (1. - pow(e_m, 2.)));
  const auto recip_d =
	1. / c_m + a_moon * e_m * cos(sm_rad - pm_rad)
	+ a_moon * pow(e_m, 2.) * cos(2. * (sm_rad - pm_rad))
	+ (15. / 8.) * a_moon * m_s2m * e_m * cos(sm_rad - 2. * hs_rad + pm_rad)
	+ a_moon * pow(m_s2m, 2.) * cos(2. * (sm_rad - hs_rad));
  return 1. / recip_d;
}

floating_t longman::distance_center_sun_earth() const noexcept { // D
  using namespace constant;
  const auto a_sun = 1. / (c_s * (1. - pow(es, 2.)));
  const auto recip_D = ((1. / c_s) + a_sun * es * cos(hs_rad - ps_rad));
  return static_cast<floating_t>(1) / recip_D;
}

floating_t longman::inclination_of_moon() const { // Im_rad
  using namespace constant;
  const auto inc_moon_rad =
	acos(cos(omega) * cos(i) - sin(omega) * sin(i) * cos(N_rad));
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

floating_t longman::longitude_celestial_equator() const { // nu
  const auto nu_t = asin(sin(constant::i) * sin(N_rad) / sin(Im_rad));
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

} // namespace iporoskun::longman
