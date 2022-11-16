#include <exception>
#include "iporoskun/longman.hpp"

namespace iporoskun::longman {

void longman::calc_longitude_and_eccentricity(floating_t T) noexcept {
  sm_rad = mean_longitude_moon(T);
  pm_rad = mean_longitude_lunar_perigee(T);
  hs_rad = mean_longitude_sun(T);
  N_rad = longitude_of_moons_ascending_node(T);
  ps_rad = mean_longitude_solar_perigee(T);
  es = eccentricity_of_earths_orbit(T);
}

auto longman::calculate_acceleration()
  -> floating_t /* meters_per_second_squared_t*/
{
  using fhours_t = std::chrono::duration<floating_t, std::ratio<3600>>;
  const auto t0 = details::from_midnight<fhours_t>(utc_time_).count();
  const auto t_rad =
	details::deg2rad(15. * (t0 - 12.) + details::rad2deg(pos_rad_cm.longitude));

  chi_m_rad = t_rad + hs_rad - nu;
  chi_s_rad = t_rad + hs_rad;

  const auto cos_alpha =
	cos(N_rad) * cos(nu) + sin(N_rad) * sin(nu) * cos(omega);
  const auto sin_alpha = (sin(omega) * sin(N_rad)) / sin(Im_rad);
  const auto alpha = 2 * atan(sin_alpha / (1. + cos_alpha));

  const auto xi = N_rad - alpha;
  const auto sigma = sm_rad - xi;

  const auto L_m_rad =
	sigma + 2. * e_m * sin(sm_rad - pm_rad)
	+ 5. / 4. * pow(e_m, 2.) * sin(2. * (sm_rad - pm_rad))
	+ 15. / 4. * m_s2m * e_m * sin(sm_rad - 2. * hs_rad + pm_rad)
	+ 11. / 8. * pow(m_s2m, 2.) * sin(2. * (sm_rad - hs_rad));

  const auto L_sun_rad = hs_rad + 2. * es * sin(hs_rad - ps_rad);

  // Zenith angle of the Moon
  const auto cos_Zm =
	sin(pos_rad_cm.latitude) * sin(Im_rad) * sin(L_m_rad)
	+ cos(pos_rad_cm.latitude)
		* (pow(cos(Im_rad / 2.), 2.) * cos(L_m_rad - chi_m_rad) + pow(sin(Im_rad / 2.), 2.) * cos(L_m_rad + chi_m_rad));

  // Zenith angle of the Sun
  const auto cos_Zs =
	sin(pos_rad_cm.latitude) * sin(omega) * sin(L_sun_rad)
	+ cos(pos_rad_cm.latitude)
		* (pow(cos(omega / 2.), 2.) * cos(L_sun_rad - chi_s_rad) + pow(sin(omega / 2.), 2.) * cos(L_sun_rad + chi_s_rad));

  // vertical component of tidal acceleration due to the Moon
  const auto gm = G * M_m * r * pow(1. / d, 3.) * (3. * pow(cos_Zm, 2.) - 1.)
				  + 1.5 * G * M_m * pow(r, 2.) * pow(1. / d, 4.)
					  * (5. * pow(cos_Zm, 3.) - 3. * cos_Zm);

  // vertical component of tidal acceleration due to the Sun :
  const auto gs = G * M_s * r * pow(1. / D, 3.) * (3. * pow(cos_Zs, 2.) - 1.);

  const auto g0_gal = (gm + gs) * beta;
  // return { g0_gal * 1000., g0_gal / 100. };
  return /*meters_per_second_squared_t*/ floating_t{ g0_gal / 100. };
}

floating_t longman::mean_longitude_moon(floating_t T) noexcept { // sm_rad
  // const auto sm_rad = details::dms2rad(270., 26., 14.72) +
  // details::deg2rad((1336. * rev_sec + 1108411.20) / 3600.)*T +
  // details::deg2rad(9.09 / 3600.)*std::pow(T, 2.) + details::deg2rad(0.0068
  // / 3600.)*std::pow(T, 3.);

  const auto sm_rad =
	details::dms2rad(270., 26., 11.72)
	+ details::deg2rad((1336. * rev_sec + 1108406.05) / 3600.) * T
	+ details::deg2rad(7.128 / 3600.) * std::pow(T, 2.)
	+ details::deg2rad(0.0072 / 3600.) * std::pow(T, 3.);

  return sm_rad;
}

floating_t
  longman::mean_longitude_lunar_perigee(floating_t T) noexcept { // pm_rad
  // const auto pm_rad = details::dms2rad(334., 19., 40.87) +
  // details::deg2rad((11. * rev_sec + 392515.94) / 3600.)*T -
  // details::deg2rad(37.24 / 3600.)*std::pow(T, 2.) - details::deg2rad(0.045
  // / 3600.)*std::pow(T, 3.);

  const auto pm_rad =
	details::dms2rad(334., 19., 46.42)
	+ details::deg2rad((11. * rev_sec + 392522.51) / 3600.) * T
	- details::deg2rad(37.15 / 3600.) * std::pow(T, 2.)
	- details::deg2rad(0.036 / 3600.) * std::pow(T, 3.);

  return pm_rad;
}

floating_t longman::mean_longitude_sun(floating_t T) noexcept { // hs_rad
  // const auto hs_rad = details::dms2rad(279., 41., 48.04) +
  // details::deg2rad(129602768.13 / 3600.)*T + details::deg2rad(1.089 /
  // 3600.)*std::pow(T, 2.);

  const auto hs_rad = details::dms2rad(279., 41., 48.05)
					  + details::deg2rad(129602768.11 / 3600.) * T
					  + details::deg2rad(1.080 / 3600.) * std::pow(T, 2.);

  return hs_rad;
}

floating_t
  longman::longitude_of_moons_ascending_node(floating_t T) noexcept { // N_rad
  // const auto N_rad = details::dms2rad(259., 10., 57.12) -
  // details::deg2rad((5. * rev_sec + 482912.63) / 3600.)*T +
  // details::deg2rad(7.58 / 3600.)*std::pow(T, 2.) + details::deg2rad(0.008 /
  // 3600.)*std::pow(T, 3.);

  const auto N_rad = details::dms2rad(259., 10., 59.81)
					 - details::deg2rad((5. * rev_sec + 482911.24) / 3600.) * T
					 + details::deg2rad(7.48 / 3600.) * std::pow(T, 2.)
					 + details::deg2rad(0.007 / 3600.) * std::pow(T, 3.);
  return N_rad;
}

floating_t
  longman::mean_longitude_solar_perigee(floating_t T) noexcept { // ps_rad
  // const auto ps_rad = details::dms2rad(281., 13., 15.00) +
  // details::deg2rad(6189.03 / 3600.)*T + details::deg2rad(1.63 /
  // 3600.)*std::pow(T, 2.) + details::deg2rad(0.012 / 3600.)*
  // std::pow(T, 3.);

  const auto ps_rad = details::dms2rad(281., 13., 14.99)
					  + details::deg2rad(6188.47 / 3600.) * T
					  + details::deg2rad(1.62 / 3600.) * std::pow(T, 2.)
					  + details::deg2rad(0.011 / 3600.) * std::pow(T, 3.);

  return ps_rad;
}

floating_t longman::eccentricity_of_earths_orbit(floating_t T) noexcept { // es
  return 0.01675104 - 0.00004180 * T - 0.000000126 * std::pow(T, 2.);
}

floating_t longman::distance_parameter(
  const longman_parameter::position_t& pos_rad_cm) noexcept { // r
  const auto C_2 =
	1. / (1. + e_crt_2 * pow(sin(pos_rad_cm.latitude.get()), 2.));
  const auto C = sqrt(C_2);
  return C * a + pos_rad_cm.height.get();
}

floating_t longman::distance_center_moon_earth() noexcept { // d
  const auto a_moon = 1. / (c_m * (1. - pow(e_m, 2.)));
  const auto recip_d =
	1. / c_m + a_moon * e_m * cos(sm_rad - pm_rad)
	+ a_moon * pow(e_m, 2.) * cos(2. * (sm_rad - pm_rad))
	+ (15. / 8.) * a_moon * m_s2m * e_m * cos(sm_rad - 2. * hs_rad + pm_rad)
	+ a_moon * pow(m_s2m, 2.) * cos(2. * (sm_rad - hs_rad));
  return 1. / recip_d;
}

floating_t longman::distance_center_sun_earth() noexcept { // D
  const auto a_sun = 1. / (c_s * (1. - pow(es, 2.)));
  const auto recip_D = ((1. / c_s) + a_sun * es * cos(hs_rad - ps_rad));
  return 1. / recip_D;
}

floating_t longman::inclination_of_moon() const { // Im_rad
  const auto inc_moon_rad =
	acos(cos(omega) * cos(i) - sin(omega) * sin(i) * cos(N_rad));
  if (
	((inc_moon_rad * 180 / std::numbers::pi) < 18)
	|| ((inc_moon_rad * 180 / std::numbers::pi) > 28)) {
	throw std::logic_error("Inclination of the Moon's orbit to the equator "
						   "is outside of the normal range.");
  }
  return inc_moon_rad;
}

floating_t longman::longitude_celestial_equator() const { // nu
  const auto nu_t = asin(sin(i) * sin(N_rad) / sin(Im_rad));
  if (
	((nu_t * 180 / std::numbers::pi) < -15)
	|| ((nu_t * 180 / std::numbers::pi) > 15)) {
	throw std::logic_error(
	  "Longitude in the celestial equator is not correct.");
  }
  return nu_t;
}

} // namespace iporoskun::longman
