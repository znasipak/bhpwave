#include "trajectory.hpp"

#define ALPHA_MIN 0.
#define ALPHA_MAX 1.
#define A_MAX 0.9999
#define OMEGA_MIN 2.e-3
//////////////////////////////
// Circular Geodesic Orbits //
//////////////////////////////

// simple signnum function for all ordered argument types
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double kerr_geo_energy_circ(double a, double r){
	double v = 1./sqrt(r);
	return (1. - 2.*pow(v, 2) + a*pow(v, 3))/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
}
double kerr_geo_momentum_circ(double a, double r){
	double v = 1./sqrt(r);
	return a*v*r*(1. - 2.*a*pow(v, 3) + pow(a, 2)*pow(v, 4))/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3))/abs(a);
}
double kerr_geo_time_frequency_circ(double a, double r){
	double v = 1./sqrt(r);
	return pow(r, 2)*(1 + a*pow(v, 3))/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
}
double kerr_geo_azimuthal_frequency_circ(double a, double r){
	double v = 1./sqrt(r);
	return pow(r, 2)*pow(v, 3)/sqrt(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3));
}

double kerr_geo_azimuthal_frequency_circ_time(double a, double r){
	double v = 1./sqrt(r);
	return pow(v, 3)/(1 + a*pow(v, 3));
}
double kerr_geo_azimuthal_frequency_circ_time(double a, double r, int sgnX){
	return kerr_geo_azimuthal_frequency_circ_time(sgnX*a, r);
}
double kerr_geo_radius_circ(double a, double Omega){
	double sgnOmega = Omega/abs(Omega);
	return pow((1. - a*Omega)/(sgnOmega*Omega), 2./3.);
}

double kerr_isco_radius(double a, int sgnX){
	double z1 = 1 + pow(1 - a*a, 1./3.)*(pow(1 - a, 1./3.) + pow(1 + a, 1./3.));
	double z2 = sqrt(3*a*a + z1*z1);

	return 3 + z2 - sgnX*sqrt((3. - z1)*(3. + z1 + 2.*z2));
}

double kerr_isco_frequency(double a){
	int sgnX = sgn(a);
	return kerr_geo_azimuthal_frequency_circ_time(abs(a), kerr_isco_radius(abs(a), sgnX), sgnX);
}

//////////////////////////////////////
// Evolution of geodesic quantities //
//////////////////////////////////////

double kerr_geo_denergy_domega_circ(double a, double om){
	return denergy_dr(a, kerr_geo_radius_circ(a, om))*dr_domega(a, om);
}
double kerr_geo_dmomentum_domega_circ(double a, double om){
	return dmomentum_dr(a, kerr_geo_radius_circ(a, om))*dr_domega(a, om);
}

double dr_domega(double a, double om){
	return -2./(3.*pow(om, 5./3.)*pow(1. - a*om, 1./3.));
}

double denergy_dr(double a, double r){
	double v = 1./sqrt(r);
	return 0.5*(pow(v, 4) - 6.*pow(v, 6) + 8.*a*pow(v, 7) - 3.*a*a*pow(v, 8))/pow(1. + v*v*(2.*a*v - 3.), 1.5);
}

double dmomentum_dr(double a, double r){
	double v = 1./sqrt(r);
	return 0.5*v*(1 + a*pow(v, 3))*(1. - 6.*pow(v,2) + 8.*a*pow(v,3) - 3.*pow(a,2)*pow(v,4))/pow(1. - 3.*pow(v, 2) + 2.*a*pow(v, 3), 1.5);
}

////////////////////////////////////////
// Quasi-circular inspiral trajectory //
////////////////////////////////////////

// double alpha_of_a_omega(const double &a, const double &omega, const double &oISCO){
// 	return pow(omega*(1. - a*omega)/oISCO/(1. - a*oISCO), 1./3.);
// }
//
// double alpha_of_a_omega(const double &a, const double &omega){
// 	return alpha_of_a_omega(a, omega, abs(kerr_isco_frequency(a)));
// }
//
// double dalpha_domega_of_a_omega(const double &a, const double &omega, const double &oISCO){
// 	return (1. - 2.*a*omega)*pow(oISCO*(1. - a*oISCO)*pow(omega*(1. - a*omega), 2), -1./3.)/3.;
// }
//
// double omega_of_a_alpha(const double &a, const double &alpha, const double &oISCO){
// 	if(a == 0){
// 		return pow(alpha, 3)*oISCO;
// 	}
// 	return 0.5*(1. - sqrt(1. - 4.*a*pow(alpha, 3)*oISCO*(1. - a*oISCO)))/a;
// }
//
// double omega_of_a_alpha(const double &a, const double &alpha){
// 	return omega_of_a_alpha(a, alpha, abs(kerr_isco_frequency(a)));
// }
//
// double spin_of_chi(const double &chi){
// 	return 1. - pow(chi, 3);
// }
//
// double chi_of_spin(const double &a){
// 	return pow(1. - a, 1./3.);
// }

double spin_of_chi_subfunc(const double & chi){
	return 1. - pow(chi, 3);
}
double chi_of_spin_subfunc(const double & a){
	return pow(1. - a, 1./3.);
}

double spin_of_chi(const double & chi){
	return 1. - pow(chi_of_spin_subfunc(A_MAX) + pow(chi, 2)*(chi_of_spin_subfunc(-A_MAX) - chi_of_spin_subfunc(A_MAX)), 3);
}

double chi_of_spin(const double & a){
	return pow((chi_of_spin_subfunc(a) - chi_of_spin_subfunc(A_MAX))/(chi_of_spin_subfunc(-A_MAX) - chi_of_spin_subfunc(A_MAX)), 0.5);
}

double alpha_of_a_omega(const double &, const double & omega, const double & oISCO){
	if(abs(oISCO - omega) < 1.e-13){return 0.;}
	return pow(abs(pow(oISCO, 1./3.) - pow(omega, 1./3.))/(pow(oISCO, 1./3.) - pow(OMEGA_MIN, 1./3.)), 0.5);
}

double alpha_of_a_omega(const double & a, const double & omega){
	return alpha_of_a_omega(a, omega, abs(kerr_isco_frequency(a)));
}

double omega_of_a_alpha(const double &, const double & alpha, const double & oISCO){
	return pow(pow(oISCO, 1./3.) - pow(alpha, 2.)*(pow(oISCO, 1./3.) - pow(OMEGA_MIN, 1./3.)), 3.);
}

double omega_of_a_alpha(const double & a, const double & alpha){
	return omega_of_a_alpha(a, alpha, abs(kerr_isco_frequency(a)));
}

double omega_of_chi_alpha(const double & chi, const double & alpha, const double & oISCO){
	return omega_of_a_alpha(spin_of_chi(chi), alpha, oISCO);
}

double omega_of_chi_alpha(const double & chi, const double & alpha){
	return omega_of_chi_alpha(chi, alpha, abs(kerr_isco_frequency(spin_of_chi(chi))));
}

double domega_dalpha_of_a_omega(const double &, const double &omega, const double &oISCO){
	if(abs(oISCO - omega) < 1.e-13){return 0.;}
	return -6.*pow((pow(oISCO, 1./3.) - pow(OMEGA_MIN, 1./3.))*(pow(oISCO, 1./3.) - pow(omega, 1./3.)), 0.5)*pow(omega, 2./3.);
}

double dalpha_domega_of_a_omega(const double &a, const double &omega, const double &oISCO){
	return 1./domega_dalpha_of_a_omega(a, omega, oISCO);
}

double max_orbital_frequency(const double &a){
  return omega_of_a_alpha(a, ALPHA_MAX);
}

double min_orbital_frequency(const double &a){
  return omega_of_a_alpha(a, ALPHA_MIN);
}

double max_orbital_radius(const double &a){
  return kerr_geo_radius_circ(a, min_orbital_frequency(a));
}

double min_orbital_radius(const double &a){
  return kerr_geo_radius_circ(a, max_orbital_frequency(a));
}

double newtonian_energy_flux(const double &omega){
  return 32./5.*pow(abs(omega), 10./3.);
}

// Trajectory Class

// Trajectory read_trajectory_data(std::string filename){
// 	double chi, alpha, t, phi, flux;
// 	Vector chiVec, alphaVec, tVec, phiVec, fluxVec;
// 	std::istringstream lin;
// 	std::ifstream inFile(filename);
// 	for (std::string line; std::getline(inFile, line); ) {
// 	    lin.clear();
// 	    lin.str(line);
// 	    if(lin >> chi >> alpha >> t >> phi >> flux){
//         chiVec.push_back(chi);
//         alphaVec.push_back(alpha);
//         tVec.push_back(t);
//         phiVec.push_back(phi);
//         fluxVec.push_back(flux);
// 	    }
// 	}
//   if(alphaVec[0] > alphaVec[1]){
//     std::reverse(chiVec.begin(), chiVec.end());
//     std::reverse(alphaVec.begin(), alphaVec.end());
//     std::reverse(tVec.begin(), tVec.end());
//     std::reverse(phiVec.begin(), phiVec.end());
// 		std::reverse(fluxVec.begin(), fluxVec.end());
//   }
//
//   Trajectory traj = {
//     .chi = chiVec,
//     .alpha = alphaVec,
//     .t = tVec,
//     .phi = phiVec,
//     .flux = fluxVec
//   };
//
// 	return traj;
// }

Data read_data(const std::string& filename){
	double x, y, z;
	Vector xVec, yVec, zVec;
	std::istringstream lin;
	std::ifstream inFile(filename);
	for (std::string line; std::getline(inFile, line); ) {
	    lin.clear();
	    lin.str(line);
	    if(lin >> x >> y >> z){
					xVec.push_back(x);
					yVec.push_back(y);
					zVec.push_back(z);
	    }
	}
	int n = xVec.size();
	EigenArray xArray(n), yArray(n), zArray(n);
	for(int i = 0; i < n; i++){
		xArray[i] = xVec[i];
		yArray[i] = yVec[i];
		zArray[i] = zVec[i];
	}

	Data data = {
		.x = xArray,
		.y = yArray,
		.z = zArray
	};

	return data;
}

TrajectoryData::TrajectoryData(const EigenArray &chi, const EigenArray &alpha, const EigenArray &t, const EigenArray &phi, const EigenArray & flux, const EigenArray &beta, const EigenArray &omega, const EigenArray &alphaOfT, const EigenArray &tMax): chi(chi), alpha(alpha), t(t), phi(phi), flux(flux), beta(beta), omega(omega), alphaOfT(alphaOfT), tMax(tMax) {}

TrajectoryData read_trajectory_data(std::string filename){
	int chiSample, alphaSample;
	double chi, alpha, t, phi, flux, beta, omega;
	// Vector chiVec, alphaVec, tVec, phiVec, fluxVec, betaVec, omegaVec, alphaOfTVec;

	std::istringstream lin;
	std::ifstream inFile(filename);
	std::string line;
	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> chiSample >> alphaSample;
	int n = chiSample*alphaSample;
	EigenArray chiA(n), alphaA(n), tA(n), phiA(n), fluxA(n), betaA(n), omegaA(n), alphaOfTA(n);
	int i = 0;

	for(std::string line; std::getline(inFile, line);){
		lin.clear();
		lin.str(line);
		if(lin >> chi >> alpha >> flux >> t >> phi >> beta >> omega){
			chiA[i] = chi;
			alphaA[i] = alpha;
			tA[i] = sqrt(log1p(t));
			phiA[i] = sqrt(log1p(phi));
			fluxA[i] = flux;
			betaA[i] = beta;
			omegaA[i] = omega;
			alphaOfTA[i] = alpha_of_a_omega(spin_of_chi(chi), omega);
			i++;
		}
	}
	EigenArray alphaAReduce(alphaSample);
	for(int i = 0; i < alphaSample; i++){
		alphaAReduce[i] = alphaA[i];
	}
	EigenArray chiAReduce(chiSample);
	for(int i = 0; i < chiSample; i++){
		chiAReduce[i] = chiA[i*alphaSample];
	}
	EigenArray tMax(chiSample);
	for(int i = 0; i < chiSample; i++){
		tMax[i] = betaA[(i + 1)*alphaSample - 1];
	}
	EigenArray betaAReduce(alphaSample);
	for(int i = 0; i < alphaSample; i++){
		betaAReduce[i] = betaA[i]/tMax[0];
	}

	TrajectoryData traj(chiAReduce, alphaAReduce, tA, phiA, fluxA, betaAReduce, omegaA, alphaOfTA, tMax);

	return traj;
}

CubicInterpolator read_time_spline(double){
  TrajectoryData traj = read_trajectory_data();
  return CubicInterpolator(traj.alpha, traj.t);
}

CubicInterpolator read_phase_spline(double){
  TrajectoryData traj = read_trajectory_data();
  return CubicInterpolator(traj.alpha, traj.phi);
}

CubicInterpolator read_flux_spline(double){
  TrajectoryData traj = read_trajectory_data();
  return CubicInterpolator(traj.alpha, traj.flux);
}

double slow_time_of_omega(CubicInterpolator &t, double a, double omega){
  return -t.evaluate(alpha_of_a_omega(a, omega));
}

double slow_phase_of_omega(CubicInterpolator &phi, double a, double omega){
  return -phi.evaluate(alpha_of_a_omega(a, omega));
}

double normalized_flux_of_omega(CubicInterpolator &Edot, double a, double omega){
  return Edot.evaluate(alpha_of_a_omega(a, omega))*newtonian_energy_flux(omega);
}

double normalized_omega_time_derivative_of_omega(CubicInterpolator &Edot, double a, double omega){
  return -(1./kerr_geo_denergy_domega_circ(a, omega) + omega/kerr_geo_dmomentum_domega_circ(a, omega))*Edot.evaluate(alpha_of_a_omega(a, omega))*newtonian_energy_flux(omega);
}

// TrajectorySpline Class

TrajectorySpline::TrajectorySpline(TrajectoryData traj):
_spin(spin_of_chi(traj.chi[0])), _omega_isco(0.), _time_spline(traj.alpha, traj.t), _phase_spline(traj.alpha, traj.phi), _flux_spline(traj.alpha, traj.flux), _frequency_spline(traj.t, traj.alpha) {_omega_isco = kerr_isco_frequency(_spin);}
TrajectorySpline::TrajectorySpline(const EigenArray & chi, const EigenArray & alpha, const EigenArray & t, const EigenArray & phi, const EigenArray & flux):
_spin(spin_of_chi(chi[0])), _omega_isco(0.), _time_spline(alpha, t), _phase_spline(alpha, phi), _flux_spline(alpha, flux), _frequency_spline(t, alpha) {_omega_isco = kerr_isco_frequency(_spin);}
TrajectorySpline::~TrajectorySpline(){}

double TrajectorySpline::time(double alpha){
	return -expm1(pow(_time_spline.evaluate(alpha), 2));
}

double TrajectorySpline::phase(double alpha){
  	return -expm1(pow(_phase_spline.evaluate(alpha), 2));
}

double TrajectorySpline::flux(double alpha){
  	return _flux_spline.evaluate(alpha)*newtonian_energy_flux(omega_of_a_alpha(_spin, alpha, _omega_isco));
}

double TrajectorySpline::flux_norm(double alpha){
  	return _flux_spline.evaluate(alpha);
}

double TrajectorySpline::orbital_alpha(double t){
  	return _frequency_spline.evaluate(sqrt(log1p(-t)));
}

double TrajectorySpline::orbital_alpha_derivative(double t){
	double sqrtlog1mt = sqrt(log1p(-t));
  	return -0.5*_frequency_spline.derivative(sqrtlog1mt)/(1. - t)/sqrtlog1mt;
}

double TrajectorySpline::orbital_frequency_time_derivative(double alpha){
	double omega = omega_of_a_alpha(_spin, alpha, _omega_isco);
  	return -(1./kerr_geo_denergy_domega_circ(_spin, omega) + omega/(kerr_geo_dmomentum_domega_circ(_spin, omega)))*_flux_spline.evaluate(alpha)*newtonian_energy_flux(omega);
}

double TrajectorySpline::orbital_frequency_time_derivative(double alpha, double omega){
  	return -(1./kerr_geo_denergy_domega_circ(_spin, omega) + omega/(kerr_geo_dmomentum_domega_circ(_spin, omega)))*_flux_spline.evaluate(alpha)*newtonian_energy_flux(omega);
}

double TrajectorySpline::time_of_omega(double omega){
  	return time(alpha_of_a_omega(_spin, omega, _omega_isco));
}

double TrajectorySpline::time_of_omega_derivative(double omega){
	double f = _time_spline.evaluate(alpha_of_a_omega(_spin, omega, _omega_isco));
  	return -2.*f*exp(pow(f, 2))*_time_spline.derivative(alpha_of_a_omega(_spin, omega, _omega_isco))*dalpha_domega_of_a_omega(_spin, omega, _omega_isco);
}

double TrajectorySpline::phase_of_omega(double omega){
  return phase(alpha_of_a_omega(_spin, omega, _omega_isco));
}

double TrajectorySpline::phase_of_omega_derivative(double omega){
	double f = _phase_spline.evaluate(alpha_of_a_omega(_spin, omega, _omega_isco));
  	return -2.*f*exp(pow(f, 2))*_phase_spline.derivative(alpha_of_a_omega(_spin, omega, _omega_isco))*dalpha_domega_of_a_omega(_spin, omega, _omega_isco);
}

double TrajectorySpline::flux_of_omega(double omega){
  	return _flux_spline.evaluate(alpha_of_a_omega(_spin, omega, _omega_isco))*newtonian_energy_flux(omega);
}

double TrajectorySpline::orbital_frequency(double t){
  	return omega_of_a_alpha(_spin, orbital_alpha(t), _omega_isco);
}

double TrajectorySpline::orbital_frequency_derivative(double t){
  	return orbital_alpha_derivative(t)/dalpha_domega_of_a_omega(_spin, orbital_frequency(t), _omega_isco);
}

double TrajectorySpline::orbital_frequency_time_derivative_of_omega(double omega){
  	return -(1./kerr_geo_denergy_domega_circ(_spin, omega) + omega/kerr_geo_dmomentum_domega_circ(_spin, omega))*_flux_spline.evaluate(alpha_of_a_omega(_spin, omega, _omega_isco))*newtonian_energy_flux(omega);
}

CubicInterpolator TrajectorySpline::get_phase_spline(){
	CubicInterpolator phase_copy(_phase_spline);
	return phase_copy;
}

double TrajectorySpline::get_spin(){
	return _spin;
}
double TrajectorySpline::get_orbital_frequency_isco(){
	return _omega_isco;
}

double TrajectorySpline::get_max_orbital_frequency(double a){
	return omega_of_a_alpha(a, ALPHA_MIN);
}

double TrajectorySpline::get_min_orbital_frequency(double a){
	return omega_of_a_alpha(a, ALPHA_MAX);
}

// TrajectorySpline2D class

TrajectorySpline2D::TrajectorySpline2D(TrajectoryData traj):
 _time_spline(traj.chi, traj.alpha, traj.t), _phase_spline(traj.chi, traj.alpha, traj.phi), _flux_spline(traj.chi, traj.alpha, traj.flux), _alpha_spline(traj.chi, traj.beta, traj.alphaOfT), _frequency_spline(traj.chi, traj.beta, traj.omega), _max_time_spline(traj.chi, traj.tMax) {}
TrajectorySpline2D::TrajectorySpline2D(const EigenArray & chi, const EigenArray & alpha, const EigenArray & beta, const EigenArray & t, const EigenArray & phi, const EigenArray & flux, const EigenArray & omega, const EigenArray & alphaOfT, const EigenArray & tMax):
  _time_spline(chi, alpha, t), _phase_spline(chi, alpha, phi), _flux_spline(chi, alpha, flux), _alpha_spline(chi, beta, alphaOfT), _frequency_spline(chi, beta, omega), _max_time_spline(chi, tMax) {}
TrajectorySpline2D::~TrajectorySpline2D(){}

double TrajectorySpline2D::time(double chi, double alpha){
  	return -expm1(pow(_time_spline.evaluate(chi, alpha), 2));
}

double TrajectorySpline2D::phase(double chi, double alpha){
  	return -expm1(pow(_phase_spline.evaluate(chi, alpha), 2));
}

double TrajectorySpline2D::flux(double chi, double alpha){
	double a = spin_of_chi(chi);
  	return _flux_spline.evaluate(chi, alpha)*newtonian_energy_flux(omega_of_a_alpha(a, alpha));
}

double TrajectorySpline2D::flux_norm(double chi, double alpha){
	double a = spin_of_chi(chi);
  	return _flux_spline.evaluate(chi, alpha);
}

double TrajectorySpline2D::orbital_alpha(double chi, double t){
  	return _alpha_spline.evaluate(chi, sqrt(log1p(-t))/_max_time_spline.evaluate(chi));
}

double TrajectorySpline2D::orbital_alpha_derivative(double chi, double t){
	double sqrtlog1mt = sqrt(log1p(-t));
	double maxt = _max_time_spline.evaluate(chi);
  	return -0.5*_alpha_spline.derivative_y(chi, sqrtlog1mt/maxt)/(1. - t)/sqrtlog1mt/maxt;
}

double TrajectorySpline2D::orbital_frequency_time_derivative(double chi, double alpha){
	double a = spin_of_chi(chi);
	double omega = omega_of_a_alpha(a, alpha);
  	return -(1./kerr_geo_denergy_domega_circ(a, omega) + omega/(kerr_geo_dmomentum_domega_circ(a, omega)))*_flux_spline.evaluate(chi, alpha)*newtonian_energy_flux(omega);
}

double TrajectorySpline2D::time_of_a_omega(double a, double omega){
  	return time(chi_of_spin(a), alpha_of_a_omega(a, omega));
}

double TrajectorySpline2D::time_of_a_omega_derivative(double a, double omega){
	double oISCO = abs(kerr_isco_frequency(a));
	double chi = chi_of_spin(a);
	double alpha = alpha_of_a_omega(a, omega, oISCO);
	double f = _time_spline.evaluate(chi, alpha);
  	return -2.*f*exp(pow(f, 2))*_time_spline.derivative_y(chi, alpha)*dalpha_domega_of_a_omega(a, omega, oISCO);
}

double TrajectorySpline2D::phase_of_a_omega(double a, double omega){
	double oISCO = abs(kerr_isco_frequency(a));
  	return phase(chi_of_spin(a), alpha_of_a_omega(a, omega));
}

double TrajectorySpline2D::phase_of_a_omega_derivative(double a, double omega){
	double oISCO = abs(kerr_isco_frequency(a));
	double chi = chi_of_spin(a);
	double alpha = alpha_of_a_omega(a, omega, oISCO);
	double f = _phase_spline.evaluate(chi, alpha);
  	return -2.*f*exp(pow(f, 2))*_phase_spline.derivative_y(chi, alpha)*dalpha_domega_of_a_omega(a, omega, oISCO);
}

double TrajectorySpline2D::flux_of_a_omega(double a, double omega){
  	return _flux_spline.evaluate(chi_of_spin(a), alpha_of_a_omega(a, omega))*newtonian_energy_flux(omega);
}

double TrajectorySpline2D::orbital_frequency(double a, double t){
  	return _frequency_spline.evaluate(chi_of_spin(a), sqrt(log(1. - t))/_max_time_spline.evaluate(chi_of_spin(a)));
}

double TrajectorySpline2D::orbital_frequency_derivative(double a, double t){
	double sqrtlog1mt = sqrt(log1p(-t));
	double chi = chi_of_spin(a);
	double maxt = _max_time_spline.evaluate(chi);
  	return -0.5*_frequency_spline.derivative_y(chi, sqrtlog1mt/maxt)/(1. - t)/sqrtlog1mt/maxt;
}

double TrajectorySpline2D::orbital_frequency_time_derivative_of_a_omega(double a, double omega){
  return -(1./kerr_geo_denergy_domega_circ(a, omega) + omega/kerr_geo_dmomentum_domega_circ(a, omega))*_flux_spline.evaluate(chi_of_spin(a), alpha_of_a_omega(a, omega))*newtonian_energy_flux(omega);
}

double TrajectorySpline2D::orbital_frequency_isco(double chi){
	return kerr_isco_frequency(spin_of_chi(chi));
}

double TrajectorySpline2D::orbital_frequency_isco_of_a(double a){
	return kerr_isco_frequency(a);
}

double TrajectorySpline2D::max_orbital_frequency(double a){
	return omega_of_a_alpha(a, ALPHA_MIN);
}

double TrajectorySpline2D::min_orbital_frequency(double a){
	return omega_of_a_alpha(a, ALPHA_MAX);
}

double TrajectorySpline2D::max_orbital_radius(double a){
	return kerr_geo_radius_circ(a, omega_of_a_alpha(a, ALPHA_MAX));
}

double TrajectorySpline2D::max_time_before_merger(double a){
	return -expm1(pow(_max_time_spline.evaluate(chi_of_spin(a)), 2));
}

BicubicInterpolator TrajectorySpline2D::get_phase_spline(){
  BicubicInterpolator phase_copy(_phase_spline);
  return phase_copy;
}
