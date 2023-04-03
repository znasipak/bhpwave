#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include "spline.hpp"
#include "omp.h"

typedef struct DataStruct{
	Vector x;
	Vector y;
	Vector z;
} Data;

class TrajectoryData{
public:
	TrajectoryData() {}
	TrajectoryData(const Vector &chi, const Vector &alpha, const Vector &t, const Vector &phi, const Vector & flux, const Vector &beta, const Vector &omega, const Vector &alphaOfT, const Vector &tMax);
	Vector chi;
	Vector alpha;
	Vector t;
  	Vector phi;
  	Vector flux;
	Vector beta;
	Vector omega;
	Vector alphaOfT;
	Vector tMax;
};

// Trajectory read_trajectory_data(std::string filename);
Data read_data(const std::string& filename);
TrajectoryData read_trajectory_data(std::string filename="data/trajectory.txt");

class TrajectorySpline{
public:
	TrajectorySpline(TrajectoryData traj);
	TrajectorySpline(const double & chi, const Vector & alpha, const Vector & t, const Vector & phi, const Vector & flux);
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

  	Spline get_phase_spline();
	double get_spin();
	double get_orbital_frequency_isco();
	double get_min_orbital_frequency(double a);
	double get_max_orbital_frequency(double a);

private:
  	double _spin;
	double _omega_isco;
  	Spline _time_spline;
  	Spline _phase_spline;
  	Spline _flux_spline;
  	Spline _frequency_spline;
};

class SmallTrajectorySpline2D{
public:
	SmallTrajectorySpline2D(std::string filename="data/trajectory.txt");
	SmallTrajectorySpline2D(TrajectoryData traj);
  	SmallTrajectorySpline2D(const Vector & chi, const Vector & alpha, const Vector & flux);
  	~SmallTrajectorySpline2D();

  	double flux(double chi, double alpha);
  	double flux_of_a_omega(double a, double omega);

private:
  	Spline2D _flux_spline;
};

class TrajectorySpline2D{
public:
	TrajectorySpline2D(std::string filename="data/trajectory.txt");
	TrajectorySpline2D(TrajectoryData traj);
  	TrajectorySpline2D(const Vector & chi, const Vector & alpha, const Vector & beta, const Vector & t, const Vector & phi, const Vector & flux, const Vector & omega, const Vector & alphaOfT, const Vector & tMax);
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

	void flux_of_a_omega(double flux[], const double a[], const double omega[], int n, int num_threads=0);

private:
  	Spline2D _time_spline;
  	Spline2D _phase_spline;
  	Spline2D _flux_spline;
  	Spline2D _alpha_spline;
	Spline2D _frequency_spline;
	Spline _max_time_spline;
};

class InspiralContainer{
public:
	InspiralContainer(int inspiralSteps);
	void setInspiralInitialConditions(double a, double massratio, double r0, double dt);
	void setTimeStep(int i, double alpha, double phase);
	
	const Vector& getAlpha() const;
	const Vector& getPhase() const;

	Vector& getAlphaNonConstRef();
	Vector& getPhaseNonConstRef();

	double getAlpha(int i);
	double getPhase(int i);
	double getTime(int i);
	double getFrequency(int i);
	double getRadius(int i);

	double getSpin();
	double getMassRatio();
	double getInitialRadius();
	double getInitialFrequency();
	int getSize();

private:
	double _a;
	double _massratio;
	double _r0;
	double _dt;
	Vector _alpha;
	Vector _phase;
};

class InspiralGenerator{
public:
	InspiralGenerator(TrajectorySpline2D &traj, int num_threads=0);
	InspiralContainer computeInspiral(double a, double massratio, double r0, double dt, double T, int num_threads=0);
	void computeInitialConditions(double &chi, double &omega_i, double &alpha_i, double &t_i, double a, double massratio, double r0, double &T);
	double computeTimeToMerger(double a, double massratio, double r0);
	int computeTimeStepNumber(double a, double massratio, double r0, double dt, double T);
	int computeTimeStepNumber(double dt, double T);

	void computeInspiral(InspiralContainer &inspiral, double chi, double omega_i, double alpha_i, double t_i, double massratio, double dt, int num_threads=0);

protected:
	TrajectorySpline2D& _traj;
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

Spline read_time_spline(double a);
Spline read_phase_spline(double a);
Spline read_flux_spline(double a);

double slow_time_of_omega(Spline &t, double a, double omega);
double slow_phase_of_omega(Spline &phi, double a, double omega);
double normalized_flux_of_omega(Spline &Edot, double a, double omega);
double normalized_omega_time_derivative_of_omega(Spline &Edot, double a, double omega);

#endif
