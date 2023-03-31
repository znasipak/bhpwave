#ifndef WAVEFORM_HPP
#define WAVEFORM_HPP

#include "harmonics.hpp"
#include "swsh.hpp"
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_trig.h>
#include <chrono>
#include <algorithm>
#include "omp.h"

#define G_const 6.67430e-11 // m^3/kg/s^2
#define GM_const 1.32712440041279419e+20 // m^3/s^2
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

class WaveformContainer{
public:
  WaveformContainer(int timeSteps);
  WaveformContainer(double *plus_ptr, double *cross_ptr, int timeSteps);
  ~WaveformContainer();
  void setTimeStep(int i, double plus, double cross);
  void addTimeStep(int i, double plus, double cross);
  void multiplyTimeStep(int i, double plus, double cross);

  double* getPlusPointer();
  double* getCrossPointer();

  double getPlus(int i);
  double getCross(int i);

  int getSize();

protected:
  double *_plus;
  double *_cross;
  int _size;
  int _owner_flag;
};

class WaveformHarmonicOptions{
public:
  WaveformHarmonicOptions(): rescale(1.), num_threads(omp_get_max_threads()), pad_output(0) {}
  WaveformHarmonicOptions(double rescale, int num, int pad_output): rescale(rescale), num_threads(num), pad_output(pad_output) {}
  
  Complex rescale;
  int num_threads;
  int pad_output;
};

class WaveformHarmonicGenerator{
public:
  WaveformHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts = HarmonicOptions(), WaveformHarmonicOptions wOpts = WaveformHarmonicOptions());

  WaveformContainer computeWaveformHarmonic(int l, int m, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts);
  WaveformContainer computeWaveformHarmonics(int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts);

  void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi);
  void computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi);
  void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi);
  void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi);

  void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts);
  void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, HarmonicOptions hOpts);
  void computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts);
  void computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts);
  void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts);
  void computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts);
  void computeWaveformHarmonics(WaveformContainer &h, int m, HarmonicSpline2D *Alm, double plusY, double crossY, InspiralContainer &inspiral, double phi, WaveformHarmonicOptions opts);

  HarmonicSelector& getModeSelector();
  HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta);
  HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts);

  WaveformHarmonicOptions getWaveformHarmonicOptions();
  HarmonicOptions getHarmonicOptions();  

protected:
  HarmonicAmplitudes& _Alm;
  HarmonicSelector _mode_selector;
  WaveformHarmonicOptions _opts;
};

class WaveformGenerator: public WaveformHarmonicGenerator{
public:
  WaveformGenerator(TrajectorySpline2D &traj, HarmonicAmplitudes &harm, HarmonicOptions hOpts = HarmonicOptions(), WaveformHarmonicOptions wOpts = WaveformHarmonicOptions());

  double convertTime(double t, double M);
  int computeTimeStepNumber(double dt, double T);
  int computeTimeStepNumber(double M, double mu, double a, double r0, double dt, double T);
  
  void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T);
  void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, WaveformHarmonicOptions wOpts);
  void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts);
  void computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts);
  
  HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T);
  HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions opts);

  void computeWaveformSourceFrame(WaveformContainer &h, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T);
private:
  InspiralGenerator _inspiralGen;
};

void sourceAngles(double &theta, double &phi, double qS, double phiS, double qK, double phiK);
Complex polarization(double qS, double phiS, double qK, double phiK);

double solar_mass_to_seconds(double mass);
double seconds_to_solar_mass(double seconds);
double solar_mass_to_meters(double mass);
double solar_mass_to_parsecs(double mass);
double parsecs_to_solar_mass(double pc);
double seconds_to_years(double seconds);
double years_to_seconds(double years);


// Old functions

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



#endif