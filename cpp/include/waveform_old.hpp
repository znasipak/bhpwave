#ifndef WAVEFORM_HPP
#define WAVEFORM_HPP

#include "harmonics.hpp"
#include "swsh.hpp"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_trig.h>
#include <chrono>
#include "omp.h"

#define G_const 6.67430e-11 // m^3/kg/s^2
#define GM_const 1.32712440041279419e+20
#define c_const 299792458. // in m/s
#define Modot_const 1.98841e+30 // in kg
#define pc_const 3.0856775814913674e+16 // in m
#define kpc_const 1e3*pc
#define Mpc_const 1e3*kpc
#define Gpc_const 1e3*Mpc
#define yr_const 31558149.763545603 // in sec (sidereal year)

typedef std::vector<float> FloatVector;
typedef std::vector<Complex> ComplexVector;
typedef std::vector<int> List;

struct GWStrain{
  GWStrain(int N);
  FloatVector plus;
  FloatVector cross;
};

struct GWStrainFourier{
  GWStrainFourier(int N);
  FloatVector plusR;
  FloatVector plusI;
  FloatVector crossR;
  FloatVector crossI;
};

double solar_mass_to_seconds(double mass);
double seconds_to_solar_mass(double seconds);
double solar_mass_to_meters(double mass);
double solar_mass_to_parsecs(double mass);
double parsecs_to_solar_mass(double pc);
double seconds_to_years(double seconds);
double years_to_seconds(double years);

double relative_harmonic_power(int l, int m, double chi, double alphaMin, double deltaAlpha, Vector &alphaDot);
void harmonic_selection(List &maxmodes, double chi, double alphaInitial, double alphaFinal, double epsilon);
double frequency_spacing(double mass1, double mass2, double spin, double r0, double &duration_yr);
double time_spacing(double mass1, double mass2, double spin, double r0, double duration_yr);
double scale_strain_amplitude(double mass1, double distance);
double scale_fourier_amplitude(double mass1, double mass2, double distance);

////////////////////////////
// Time Domain Trajectory //
////////////////////////////

double time_to_inspiral(double mass1, double mass2, double spin, double r0);
int inspiral_time_steps(double mass1, double mass2, double spin, double r0, double duration_yr, double deltaT_sec);
void generate_trajectory_td(double mass1, double mass2, double spin, double r0, double phi0, double duration_yr, double deltaT_sec = 0.);
void output_trajectory_td(double *rp, double *phase, double mass1, double mass2, double spin, double r0, double phi0, double duration_yr, double deltaT_sec = 0.);
void output_downsampled_trajectory(double *t, double *alphaOfT, double *phaseOfT, int t_steps, double theta, double phi, double mass1, double mass2, double spin, double r0, double duration_yr, double deltaT_sec);
void output_flux(double *flux, double *a, double *r, int sampleN);

///////////////////////////
// Time Domain Waveforms //
///////////////////////////

void reduced_waveform_harmonic_td(FloatVector &hp, FloatVector &hc, double &chi, Vector &alpha, List &maxmodes, double &theta, double &phi, Vector &phase, int num_threads);
void reduced_waveform_harmonic_td(FloatVector &hp, FloatVector &hc, double &chi, Vector &alpha, int &L, int &m, double &theta, double &phi, Vector &phase, int num_threads=0);
GWStrain waveform_harmonic_td(int L, int m, double theta, double phi, double mass1, double mass2, double spin, double r0, double duration, double deltaT, int num_threads=0);
GWStrain waveform_td(double theta, double phi, double mass1, double mass2, double spin, double r0, double duration, double deltaT, int num_threads=0);
void generate_waveform_harmonic_td(int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec = 0., int num_threads=0);
void generate_waveform_td(double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec = 0., int num_threads=0);
void output_waveform_td(float *hp, float *hc, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec_temp, int num_threads=0);
void output_waveform_harmonic_td(float *hp, float *hc, int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec_temp, int num_threads=0);

////////////////////////////////
// Frequency Domain Waveforms //
////////////////////////////////

void reduced_waveform_harmonic_fd(FloatVector &hp, FloatVector &hc, double chi, const Vector &freq, const int &L, const int &m, const double &t0, const double &theta, const double &phi, const double &phase0, const double &massratio, const double &omegaMin, const double &omegaMax);
void reduced_waveform_harmonic_fd(FloatVector &hpR, FloatVector &hpI, FloatVector &hcR, FloatVector &hcI, double chi, const Vector &freq, const int &l, const int &m, const double &t0, const double &theta, const double &phi, const double &phase0, const double &massratio, const double &omegaMin, const double &omegaMax);
GWStrainFourier waveform_harmonic_fd(int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double &deltaF_nHz, double theta, double phi);
GWStrainFourier waveform_fd(double mass1, double mass2, double spin, double r0, double duration, double deltaF_nHz, double theta, double phi);
void generate_waveform_harmonic_fd(int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaF_nHz = 0.);
void generate_waveform_fd(double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaF_nHz = 0.);

//////////////////////////////
// Time-Frequency Waveforms //
//////////////////////////////

void reduced_waveform_harmonic_tf(FloatVector &hp, FloatVector &hc, double chi, const Vector &freq, const int &L, const int &m, const double &t0, const double &theta, const double &phi, const double &phase0, const double &massratio, const double &omegaMin, const double &omegaMax);

#endif