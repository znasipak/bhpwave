#ifndef FOURIER_HPP
#define FOURIER_HPP

#include "waveform.hpp"

class WaveformFourierHarmonicGenerator{
public:
  WaveformFourierHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts = HarmonicOptions(), WaveformHarmonicOptions wOpts = WaveformHarmonicOptions());

  void computeWaveformFourierHarmonics(WaveformContainer &h, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, HarmonicOptions hOpts, int num_threads, double freq[], int fsamples);
  void computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples, int include_negative_m = 1);
  void computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples);

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

class WaveformFourierGenerator: public WaveformFourierHarmonicGenerator{
public:
  WaveformFourierGenerator(TrajectorySpline2D &traj, HarmonicAmplitudes &harm, HarmonicOptions hOpts = HarmonicOptions(), WaveformHarmonicOptions wOpts = WaveformHarmonicOptions());

  double convertTime(double t, double M);
  double convertFrequency(double f, double M);
  int computeFrequencyStepNumber(double df, double T);
  int computeTimeStepNumber(double df, double T);
  int computeTimeStepNumber(double M, double mu, double a, double r0, double dt, double T);
  
  void computeFourierWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts);
  void computeFourierWaveform(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts);
  
  void computeFourierWaveformSourceFrame(WaveformContainer &h, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double frequencies[], double T);
  void computeFourierWaveformSourceFrame(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double frequencies[], double T);

  HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double T);
  HarmonicModeContainer selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double T, HarmonicOptions opts);

private:
  InspiralGenerator _inspiralGen;
};

double scale_fourier_amplitude(double mass1, double mass2, double distance);

#endif
