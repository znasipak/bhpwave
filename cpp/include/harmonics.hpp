#ifndef HARMONICS_HPP
#define HARMONICS_HPP

#include <string>
#include <utility>
#include <map>
#include "trajectory.hpp"
#include "swsh.hpp"
#include "omp.h"

typedef struct HarmonicModeStruct{
	Vector chi;
	Vector alpha;
	Vector A;
	Vector Phi;
} HarmonicModeData;

HarmonicModeData read_harmonic_mode_data(int L, int m, std::string filepath_base = "data/circ_data");

class HarmonicSpline{
public:
  HarmonicSpline(double chi, const Vector &alpha, const Vector &A, const Vector &Phi);
  HarmonicSpline(double spin, Spline amplitude_spline, Spline phase_spline);
  ~HarmonicSpline();

  double amplitude(double alpha);
  double phase(double alpha);

	double amplitude_of_omega(double omega);
  double phase_of_omega(double omega);
	double phase_of_omega_derivative(double omega);

private:
  double _spin;
  Spline _amplitude_spline;
  Spline _phase_spline;
};

class HarmonicSpline2D{
public:
  HarmonicSpline2D(int L, int m, std::string filepath_base = "data/circ_data");
  HarmonicSpline2D(HarmonicModeData mode);
  HarmonicSpline2D(const Vector &chi, const Vector &alpha, const Vector &Amp, const Vector &Phi);
  ~HarmonicSpline2D();

  double amplitude(double chi, double alpha);
  double phase(double chi, double alpha);

	double amplitude_of_a_omega(double a, double omega);
  double phase_of_a_omega(double a, double omega);
	double phase_of_a_omega_derivative(double a, double omega);

private:
  Spline2D _amplitude_spline;
  Spline2D _phase_spline;
};

class HarmonicAmplitudes{
public:
  HarmonicAmplitudes(std::vector<int> &lmodes, std::vector<int> &mmodes, std::string filepath_base = "data/circ_data");
  HarmonicAmplitudes(int lmodes[], int mmodes[], int modeNum, std::string filepath_base = "data/circ_data");
  ~HarmonicAmplitudes();

  double amplitude(int l, int m, double chi, double alpha);
  double phase(int l, int m, double chi, double alpha);
  int key_check(std::pair<int, int> key);

	double amplitude_of_a_omega(int l, int m, double a, double omega);
  double phase_of_a_omega(int l, int m, double a, double omega);
	double phase_of_a_omega_derivative(int l, int m, double a, double omega);
  
  HarmonicSpline2D* getPointer(int l, int m);

private:
  int _modeNum;
  std::vector<HarmonicSpline2D*> _harmonics;
  std::map<std::pair<int,int>,int> _position_map;
};

class HarmonicModeContainer{
public:
  HarmonicModeContainer() {}
  std::vector<int> lmodes;
  std::vector<int> mmodes;
  Vector plusY;
  Vector crossY;
};

class HarmonicOptions{
public:
  HarmonicOptions(): epsilon(1.e-5), max_samples(50)  {}
  HarmonicOptions(double eps, int max): epsilon(eps), max_samples(max)  {}
  double epsilon;
  int max_samples;
};

class HarmonicSelector{
public:
  HarmonicSelector(HarmonicAmplitudes &harm, HarmonicOptions opts = HarmonicOptions());

  double modePower(int l, int m, InspiralContainer &inspiral);
  int gradeMode(int l, int m, InspiralContainer &inspiral, double power22);
  int gradeMode(int l, int m, InspiralContainer &inspiral, double power22, double &plusYlm, double &crossYlm, double theta);
  HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta);

  double modePower(int l, int m, InspiralContainer &inspiral, HarmonicOptions opts);
  int gradeMode(int l, int m, InspiralContainer &inspiral, double power22, HarmonicOptions opts);
  int gradeMode(int l, int m, InspiralContainer &inspiral, double power22, double &plusYlm, double &crossYlm, double theta, HarmonicOptions opts);
  HarmonicModeContainer selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts);

  HarmonicOptions getHarmonicOptions();

private:
  HarmonicAmplitudes& _harm;
  HarmonicOptions _opts;
};

double Yslm_plus_polarization(int l, int m, double theta);
double Yslm_cross_polarization(int l, int m, double theta);
void Yslm_plus_cross_polarization(double &plusY, double &crossY, int l, int m, double theta);

#endif