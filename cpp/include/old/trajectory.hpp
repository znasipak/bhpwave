#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include "spline.hpp"

typedef struct DataStruct{
	EigenArray x;
	EigenArray y;
	EigenArray z;
} Data;

class TrajectoryData{
public:
	TrajectoryData() {}
	TrajectoryData(const EigenArray &chi, const EigenArray &alpha, const EigenArray &t, const EigenArray &phi, const EigenArray & flux, const EigenArray &beta, const EigenArray &omega, const EigenArray &alphaOfT, const EigenArray &tMax);
	EigenArray chi;
	EigenArray alpha;
	EigenArray t;
  	EigenArray phi;
  	EigenArray flux;
	EigenArray beta;
	EigenArray omega;
	EigenArray alphaOfT;
	EigenArray tMax;
};

// Trajectory read_trajectory_data(std::string filename);
Data read_data(const std::string& filename);
TrajectoryData read_trajectory_data(std::string filename="data/trajectory.txt");

class TrajectorySpline{
public:
	TrajectorySpline(TrajectoryData traj);
	TrajectorySpline(const EigenArray & chi, const EigenArray & alpha, const EigenArray & t, const EigenArray & phi, const EigenArray & flux);
	~TrajectorySpline();

	double time(double alpha);
	double phase(double alpha);
	double flux(double alpha);
	double flux_norm(double alpha);
	double orbital_alpha(double t);
	double orbital_alpha_derivative(double t);
	double orbital_frequency_time_derivative(double alpha);
	double orbital_frequency_time_derivative(double alpha, double omega);

  	double time_of_omega(double omega);
	double time_of_omega_derivative(double omega);
  	double phase_of_omega(double omega);
	double phase_of_omega_derivative(double omega);
  	double flux_of_omega(double omega);
	double orbital_frequency(double t);
	double orbital_frequency_derivative(double t);
  	double orbital_frequency_time_derivative_of_omega(double omega);

  	CubicInterpolator get_phase_spline();
	double get_spin();
	double get_orbital_frequency_isco();
	double get_min_orbital_frequency(double a);
	double get_max_orbital_frequency(double a);

private:
  	double _spin;
	double _omega_isco;
  	CubicInterpolator _time_spline;
  	CubicInterpolator _phase_spline;
  	CubicInterpolator _flux_spline;
  	CubicInterpolator _frequency_spline;
};

class TrajectorySpline2D{
public:
	TrajectorySpline2D(TrajectoryData traj);
  	TrajectorySpline2D(const EigenArray & chi, const EigenArray & alpha, const EigenArray & beta, const EigenArray & t, const EigenArray & phi, const EigenArray & flux, const EigenArray & omega, const EigenArray & alphaOfT, const EigenArray & tMax);
  	~TrajectorySpline2D();

	double time(double chi, double alpha);
  	double phase(double chi, double alpha);
  	double flux(double chi, double alpha);
	double flux_norm(double chi, double alpha);
	double orbital_alpha(double chi, double t);
	double orbital_alpha_derivative(double chi, double t);
	double orbital_frequency_time_derivative(double chi, double alpha);

  	double time_of_a_omega(double a, double omega);
	double time_of_a_omega_derivative(double a, double omega);
  	double phase_of_a_omega(double a, double omega);
	double phase_of_a_omega_derivative(double a, double omega);
  	double flux_of_a_omega(double a, double omega);
	double orbital_frequency(double a, double t);
	double orbital_frequency_derivative(double a, double t);
  	double orbital_frequency_time_derivative_of_a_omega(double a, double omega);

	double orbital_frequency_isco(double chi);
	double orbital_frequency_isco_of_a(double a);
	double min_orbital_frequency(double a);
	double max_orbital_frequency(double a);
	double max_orbital_radius(double a);
	double max_time_before_merger(double a);

  	BicubicInterpolator get_phase_spline();

private:
  	BicubicInterpolator _time_spline;
  	BicubicInterpolator _phase_spline;
  	BicubicInterpolator _flux_spline;
  	BicubicInterpolator _alpha_spline;
	BicubicInterpolator _frequency_spline;
	CubicInterpolator _max_time_spline;
};

//////////////////////////////
// Circular Geodesic Orbits //
//////////////////////////////

double kerr_geo_energy_circ(double a, double r);
double kerr_geo_momentum_circ(double a, double r);
double kerr_geo_time_frequency_circ(double a, double r);
double kerr_geo_azimuthal_frequency_circ(double a, double r);

double kerr_geo_azimuthal_frequency_circ_time(double a, double r, int sgnX);
double kerr_geo_azimuthal_frequency_circ_time(double a, double r);
double kerr_geo_radius_circ(double a, double Omega);

double kerr_isco_radius(double a, int sgnX);
double kerr_isco_frequency(double a);

//////////////////////////////////////
// Evolution of geodesic quantities //
//////////////////////////////////////

double kerr_geo_denergy_domega_circ(double a, double om);
double kerr_geo_dmomentum_domega_circ(double a, double om);
double dr_domega(double a, double om);
double denergy_dr(double a, double r);
double dmomentum_dr(double a, double r);

////////////////////////////////////////
// Quasi-circular inspiral trajectory //
////////////////////////////////////////

double alpha_of_a_omega(const double &a, const double &omega, const double &oISCO);
double alpha_of_a_omega(const double &a, const double &omega);
double dalpha_domega_of_a_omega(const double &a, const double &omega, const double &oISCO);
double omega_of_a_alpha(const double &a, const double &alpha, const double &oISCO);
double omega_of_a_alpha(const double &a, const double &alpha);
double spin_of_chi(const double &chi);
double chi_of_spin(const double &a);

double max_orbital_frequency(const double &a);
double min_orbital_frequency(const double &a);
double max_orbital_radius(const double &a);
double min_orbital_radius(const double &a);

CubicInterpolator read_time_spline(double a);
CubicInterpolator read_phase_spline(double a);
CubicInterpolator read_flux_spline(double a);

double slow_time_of_omega(CubicInterpolator &t, double a, double omega);
double slow_phase_of_omega(CubicInterpolator &phi, double a, double omega);
double normalized_flux_of_omega(CubicInterpolator &Edot, double a, double omega);
double normalized_omega_time_derivative_of_omega(CubicInterpolator &Edot, double a, double omega);

#endif
