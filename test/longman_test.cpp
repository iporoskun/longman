#include <catch2/catch.hpp>
#include <iporoskun/longman.hpp>
#define _USE_MATH_DEFINES
#include <math.h>

// https://sbgf.org.br/revista/index.php/rbgf/article/viewFile/793/416

using namespace std::chrono;
using namespace std::chrono_literals;

namespace longman = iporoskun::longman;
using position_t = longman::longman_parameter::position_t;

auto pos = position_t{ longman::latitude(52.29610),
					   longman::longitude(10.45900), longman::height(80) };


TEST_CASE("angle transformations") {

  SECTION("deg min sec to radians") {

	constexpr auto positive = longman::details::dms_to_rad(270., 26., 11.72);
	constexpr auto negative = longman::details::dms_to_rad(334., 19., 46.42);

	CHECK(positive == Approx(4.720008893));
	CHECK(negative == Approx(5.835151628));
  }

  SECTION("degree to rad") {
	CHECK(longman::details::deg_to_rad(45) == Approx(0.785398163397));
  }

  SECTION("megree-min-sec to deg") {
	CHECK(
	  longman::degree_minute_second_to_degree(22, 44, 0)
	  == Approx(22.7333333333));
  }
}


TEST_CASE("amount of centureis from ref-date") {
  using longman::details::julian_centuries_from_ref_date;

  SECTION("dates are same ") {
	CHECK(
	  julian_centuries_from_ref_date(sys_days{ 1899y / 12 / 31d } + 12h) == 0.);
  }

  SECTION("1 century") {
	// https://www.wolframalpha.com/input/?i=31+december+1899+%2B+%28100+*+365.25+days%29+
	CHECK(
	  julian_centuries_from_ref_date(sys_days{ 2000y / 1 / 1d } + 12h) == 1.);
  }

  SECTION("matlab result 1") {
	CHECK(
	  julian_centuries_from_ref_date(sys_days{ 2010y / 10 / 31d })
	  == Approx(1.108302722640506));
  }

  SECTION("matlab result 2") {
	CHECK(
	  julian_centuries_from_ref_date(sys_days{ 2020y / 12 / 7d })
	  == Approx(1.209329416685681));
  }
}

TEST_CASE("time from midnight") {
  const auto time = sys_days{ (year_month_day{ 2010y / 9 / 30d }) };
  const auto T = time + 6h + 10min + 0s;

  using longman::details::from_midnight;
  using fhours_t = std::chrono::duration<double, std::ratio<3600>>;
  CHECK(from_midnight<fhours_t>(time).count() == 0.);
  CHECK(from_midnight<fhours_t>(T).count() == Approx(6.166666666666667));
}

TEST_CASE("creating longman object") {

  SECTION("position and utc offset-time") { longman::longman longman(pos, 1h); }
}

TEST_CASE("using longman's operator()") {

  longman::longman longman(pos, 1h);

  SECTION("") {
	auto now = std::chrono::system_clock::now();
	const auto result = longman(now);
	CHECK(result != 0);
  }
}

TEST_CASE("comparision with matlab results") {
  // Longman(-22.733,-90.50,0,31,12,1899,0,10,6,0,31,10,2010);
  const auto day = sys_days{ (year_month_day{ 2010y / 10 / 31d }) };
  const auto time = day + 6h + 10min + 0s;
  const double T = longman::details::julian_centuries_from_ref_date(time);
  const auto position =
	position_t{ longman::latitude(-22.733), longman::longitude(-90.50),
				longman::height(0) };
  longman::longman longman(position, 0h);

  const auto EPS = 1e-14;

  const auto T_matlab = 1.108302722640506;

  SECTION("time") {
	auto mins = ((T - T_matlab) * 365.25 * 24 * 60 * 100);
	const auto h = mins / 60;

	INFO("delta_time = " << h << " h ( mins = " << mins << " )");
	CHECK(T == Approx(T_matlab).epsilon(1e-10));
  }

  SECTION("longitude_and_eccentricity") {

	SECTION("T from matlab") {
	  const auto result = 9314.140709778332;
	  CHECK(longman::longman::mean_longitude_moon(T_matlab) == Approx(result));
	}

	CHECK(
	  longman::longman::mean_longitude_moon(T)
	  == Approx(9314.140709778332).epsilon(EPS));
	CHECK(
	  longman::longman::mean_longitude_lunar_perigee(T)
	  == Approx(84.544418585647875).epsilon(EPS));
	CHECK(
	  longman::longman::mean_longitude_sun(T)
	  == Approx(7.012636463119221e+02).epsilon(EPS));
	CHECK(
	  longman::longman::longitude_of_moons_ascending_node(T)
	  == Approx(-32.889490944928838).epsilon(EPS));
	CHECK(
	  longman::longman::mean_longitude_solar_perigee(T)
	  == Approx(4.941491045285590).epsilon(EPS));
	CHECK(
	  longman::longman::eccentricity_of_earths_orbit(T)
	  == Approx(0.016704558175993).epsilon(EPS));
  }

  const auto discard = longman(time);

  SECTION("distances in cm") {
	CHECK(
	  longman.distance_parameter(longman.get_pos_rad_cm())
	  == (6.375063476365359e+08));

	CHECK(
	  longman.distance_center_moon_earth()
	  == Approx(3.688271991063581e+10).margin(0.03));

	CHECK(
	  longman.distance_center_sun_earth()
	  == Approx(1.484005330928693e+13).margin(0.03));

	CHECK(
	  longman.longitude_celestial_equator()
	  == Approx(-0.217346347550123).epsilon(10 * EPS));
  }

  SECTION("angles in rad") {
	CHECK(
	  longman.inclination_of_moon()
	  == Approx(0.426733793977650).epsilon(10 * EPS));
  }

  SECTION("acceleration") {
	CHECK(longman(time) == Approx(-0.00000035704380751590).epsilon(1e-10));
  }
}
