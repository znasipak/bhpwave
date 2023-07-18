#include <iostream>
#include <vector>
#include <complex>

typedef std::vector<double> Vector;
typedef std::complex<double> Complex;

double factorial(int n);
double binomial(int n, int m);

// scalar spherical harmonics or normalized associated legendre polynomials
double spherical_harmonic(const int &l, const int &m, const double &theta);
Vector spherical_harmonic(const int &l, const int &m, const Vector &theta);
Complex spherical_harmonic(const int &l, const int &m, const double &th, const double &ph);

// spin-weighted spherical harmonics
double spin_weighted_spherical_harmonic(int s, int l, int m, double theta);
Vector spin_weighted_spherical_harmonic(int s, int l, int m, Vector theta);
Complex spin_weighted_spherical_harmonic(int s, int l, int m, double theta, double phi);
void spin_weighted_spherical_harmonic(double *yslm, int pts_num, int st, int l, int mt, double *theta);