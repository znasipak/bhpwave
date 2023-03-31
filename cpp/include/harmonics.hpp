#ifndef HARMONICS_HPP
#define HARMONICS_HPP

#include <string>
#include <utility>
#include <map>
#include "trajectory.hpp"
#include "swsh.hpp"
#include "omp.h"

typedef struct HarmonicModeStruct{
	EigenArray chi;
	EigenArray alpha;
	EigenArray A;
	EigenArray Phi;
} HarmonicModeData;

HarmonicModeData read_harmonic_mode_data(int L, int m, std::string filepath_base = "data/circ_data");

class HarmonicSpline{
public:
  HarmonicSpline(double chi, const EigenArray &alpha, const EigenArray &A, const EigenArray &Phi);
  HarmonicSpline(double spin, EigenCubicInterpolator amplitude_spline, EigenCubicInterpolator phase_spline);
  ~HarmonicSpline();

  double amplitude(double alpha);
  double phase(double alpha);

	double amplitude_of_omega(double omega);
  double phase_of_omega(double omega);
	double phase_of_omega_derivative(double omega);

private:
  double _spin;
  EigenCubicInterpolator _amplitude_spline;
  EigenCubicInterpolator _phase_spline;
};

class HarmonicSpline2D{
public:
  HarmonicSpline2D(int L, int m, std::string filepath_base = "data/circ_data");
  HarmonicSpline2D(HarmonicModeData mode);
  HarmonicSpline2D(const EigenArray &chi, const EigenArray &alpha, const EigenArray &Amp, const EigenArray &Phi);
  ~HarmonicSpline2D();

  double amplitude(double chi, double alpha);
  double phase(double chi, double alpha);

	double amplitude_of_a_omega(double a, double omega);
  double phase_of_a_omega(double a, double omega);
	double phase_of_a_omega_derivative(double a, double omega);

  // Spline2D getAmplitudeSpline();
  // Spline2D getPhaseSpline();
  // EigenCubicInterpolator getReducedAmplitudeSpline(double a);
  // EigenCubicInterpolator getReducedPhaseSpline(double a);

  // HarmonicSpline constant_spin_spline(double a);

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

// class HarmonicMode{
// public:
//   HarmonicMode(const Eigen::ArrayXi &lmodes, const Eigen::ArrayXi &mmodes);
//   void generateModeData(EigenArray alpha);
//   Eigen::ArrayXi getLModes();
//   Eigen::ArrayXi getMModes();
//   Eigen::ArrayXi getAmplitudes();
//   Eigen::ArrayXi getPhases();

// private:
//   Eigen::ArrayXi lmodes;
//   Eigen::ArrayXi mmodes;
//   Eigen::MatrixXd Aplus;
//   Eigen::MatrixXd Across;
//   Eigen::MatrixXd phases;
// };

#endif