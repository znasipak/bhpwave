#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <sys/stat.h>
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
	TrajectoryData(const Vector &chi, const Vector &alpha, const Vector &t, const Vector &phi, const Vector &chiFlux, const Vector &alphaFlux, const Vector & flux, const Vector &beta, const Vector &omega, const Vector &alphaOfT, const Vector &phiOfT, const double &tMax);
	Vector chi;
	Vector alpha;
	Vector t;
  	Vector phi;
	Vector chiFlux;
	Vector alphaFlux;
  	Vector flux;
	Vector beta;
	Vector omega;
	Vector alphaOfT;
	Vector phiOfT;
	double tMax;
};

// Trajectory read_trajectory_data(std::string filename);
Data read_data(const std::string& filename);
TrajectoryData read_trajectory_data(std::string filename="data/trajectory.txt");

class TrajectorySpline2D{
public:
	TrajectorySpline2D(std::string filename="data/trajectory.txt");
	TrajectorySpline2D(TrajectoryData traj);
  	TrajectorySpline2D(const Vector & chi, const Vector & alpha, const Vector & chiFlux, const Vector & alphaFlux, const Vector & beta, const Vector & t, const Vector & phi, const Vector & flux, const Vector & omega, const Vector & alphaOfT, const Vector & phaseOfT, const double & tMax);
  	~TrajectorySpline2D();

	// frequency domain
	double time(double chi, double alpha);
  	double phase(double chi, double alpha);
  	double flux(double chi, double alpha);
	double flux_norm(double chi, double alpha);

	double time_of_a_omega(double a, double omega);
	double time_of_a_omega_derivative(double a, double omega);
  	double phase_of_a_omega(double a, double omega);
	double phase_of_a_omega_derivative(double a, double omega);
  	double flux_of_a_omega(double a, double omega);
	double flux_of_a_omega_derivative(double a, double omega);
	double flux_of_a_derivative_omega(double a, double omega);

	double orbital_frequency_time_derivative_from_flux(double chi, double alpha);
	double orbital_frequency_time_derivative_from_flux_of_a_omega(double a, double omega);
	double time_of_a_alpha_omega_derivative(double a, double alpha);

	// time domain
	double orbital_alpha(double chi, double t);
	double orbital_alpha_derivative(double chi, double t);
	double orbital_frequency(double a, double t);
	double orbital_frequency_derivative(double a, double t);
	double phase_of_time(double chi, double t);
	double phase_of_time_derivative(double chi, double t);
	double phase_of_a_time(double a, double t);
	double phase_of_a_time_derivative(double a, double t);

	// utility
	double orbital_frequency_isco(double chi);
	double orbital_frequency_isco_of_a(double a);
	double min_orbital_frequency(double a);
	double max_orbital_frequency(double a);
	double max_orbital_radius(double a);
	double max_time_before_merger(double a);

	void flux_of_a_omega(double flux[], const double a[], const double omega[], int n, int num_threads=0);

private:
  	BicubicSpline _time_spline;
  	BicubicSpline _phase_spline;
  	BicubicSpline _flux_spline;
  	BicubicSpline _alpha_spline;
	BicubicSpline _frequency_spline;
	BicubicSpline _phase_time_spline;
	double _time_norm_parameter;
};

class InspiralContainer{
public:
	InspiralContainer(int inspiralSteps);
	void setInspiralInitialConditions(double a, double massratio, double r0, double dt);
	void setTimeStep(int i, double alpha, double phase);
	void setTimeStep(int i, double alpha, double phase, double dtdo);
	
	const Vector& getAlpha() const;
	const Vector& getPhase() const;

	Vector& getAlphaNonConstRef();
	Vector& getPhaseNonConstRef();

	double getAlpha(int i);
	double getPhase(int i);
	double getTimeOmegaDeriv(int i);
	double getTime(int i);
	double getFrequency(int i);
	double getRadius(int i);

	double getSpin();
	double getMassRatio();
	double getInitialRadius();
	double getFinalRadius();
	double getInitialFrequency();
	double getFinalFrequency();
	double getTimeSpacing();
	double getDuration();
	double getISCOFrequency();
	double getISCORadius();
	int getSize();

private:
	double _a;
	double _massratio;
	double _r0;
	double _dt;
	double _risco;
	double _oisco;
	Vector _alpha;
	Vector _phase;
};

class InspiralGenerator{
public:
	InspiralGenerator(TrajectorySpline2D &traj, int num_threads=0);
	InspiralContainer computeInspiral(double a, double massratio, double r0, double dt, double T, int num_threads = 0);
	void computeInitialConditions(double &chi, double &omega_i, double &alpha_i, double &t_i, double a, double massratio, double r0, double &T);
	double computeTimeToMerger(double a, double massratio, double r0);
	int computeTimeStepNumber(double a, double massratio, double r0, double dt, double T);
	int computeTimeStepNumber(double dt, double T);

	void computeInspiral(InspiralContainer &inspiral, double chi, double omega_i, double alpha_i, double t_i, double massratio, double dt, int num_threads=0);
	TrajectorySpline2D& getTrajectorySpline();

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
double kerr_isco_radius(double a);
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

CubicSpline read_time_spline(double a);
CubicSpline read_phase_spline(double a);
CubicSpline read_flux_spline(double a);

double slow_time_of_omega(CubicSpline &t, double a, double omega);
double slow_phase_of_omega(CubicSpline &phi, double a, double omega);
double normalized_flux_of_omega(CubicSpline &Edot, double a, double omega);
double normalized_omega_time_derivative_of_omega(CubicSpline &Edot, double a, double omega);

#endif
