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

double kerr_isco_radius(double a){
	return kerr_isco_radius(a, sgn(a));
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

double dalpha_domega_of_alpha(const double &alpha, const double &oISCO){
	double oISCOThird = pow(oISCO, 1./3.);
	double delta = oISCOThird - pow(OMEGA_MIN, 1./3.);
	return -1./(6.*alpha*delta*pow(oISCOThird - alpha*alpha*delta, 2.));
}

double beta_of_gamma(double gamma, double gammaMax){
	return gamma/gammaMax;
}

double gamma_of_time(double time){
	return pow(log1p(-time), 1./6.);
}

double beta_of_time(double time, double gammaMax){
	return gamma_of_time(time)/gammaMax;
}

double time_of_gamma(double gamma){
	return -expm1(pow(gamma, 6.));
}

double time_of_beta(double beta, double gammaMax){
	return -expm1(pow(beta*gammaMax, 6.));
}

double dbeta_dtime(double time, double gammaMax){
	return ((1./6.)/(time - 1)*pow(log1p(-time), -5./6.))/gammaMax;
}

double phase_norm_2(double phase){
	return log1p(-phase);
}

double phase_of_phase_norm_2(double phase2){
	return -expm1(phase2);
}

double dphase_dphase_norm_2(double phase2){
	return -exp(phase2);
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

double normalize_energy_flux(const double &omega){
  return pow(abs(omega), 10./3.);
}

double normalize_time(const double &omega, const double &oISCO){
  return pow(abs(omega), -8./3.) - pow(abs(oISCO), -8./3.) + 1.e-6;
}

double normalize_phase(const double &omega, const double &oISCO){
  return  pow(abs(omega), -5./3.) - pow(abs(oISCO), -5./3.) + 1.e-6;
}

double normalize_time_domega(const double &omega){
  return (-8./3.)*pow(abs(omega), -11./3.);
}

double normalize_phase_domega(const double &omega){
  return  (-5./3.)*pow(abs(omega), -8./3.);
}

// Trajectory Class

InspiralContainer::InspiralContainer(int inspiralSteps): _alpha(inspiralSteps), _phase(inspiralSteps) {}
void InspiralContainer::setInspiralInitialConditions(double a, double massratio, double r0, double dt){
  _a = a;
  _massratio = massratio;
  _r0 = r0;
  _dt = dt;
  _risco = kerr_isco_radius(_a);
  _oisco = kerr_isco_frequency(_a);
}
void InspiralContainer::setTimeStep(int i, double alpha, double phase){
  _alpha[i] = alpha;
  _phase[i] = phase;
}
void InspiralContainer::setTimeStep(int i, double alpha, double phase, double dtdo){
  _alpha[i] = alpha;
  _phase[i] = phase;
}

const Vector& InspiralContainer::getAlpha() const{
  return _alpha;
}
const Vector& InspiralContainer::getPhase() const{
  return _phase;
}

Vector& InspiralContainer::getAlphaNonConstRef(){
  return _alpha;
}
Vector& InspiralContainer::getPhaseNonConstRef(){
  return _phase;
}

double InspiralContainer::getAlpha(int i){
  return _alpha[i];
}
double InspiralContainer::getPhase(int i){
  return _phase[i];
}

double InspiralContainer::getTime(int i){
	return i*_dt;
}
double InspiralContainer::getFrequency(int i){
	return omega_of_a_alpha(_a, _alpha[i]);
}
double InspiralContainer::getRadius(int i){
	return kerr_geo_radius_circ(_a, omega_of_a_alpha(_a, _alpha[i]));
}

double InspiralContainer::getSpin(){
  return _a;
}
double InspiralContainer::getISCORadius(){
  return _risco;
}
double InspiralContainer::getISCOFrequency(){
  return _oisco;
}

double InspiralContainer::getMassRatio(){
  return _massratio;
}

double InspiralContainer::getInitialRadius(){
  return _r0;
}

double InspiralContainer::getInitialFrequency(){
  return kerr_geo_azimuthal_frequency_circ_time(_a, _r0);
}

double InspiralContainer::getFinalFrequency(){
  return omega_of_a_alpha(_a, _alpha[getSize() - 1], _oisco);
}

double InspiralContainer::getFinalRadius(){
  return kerr_geo_radius_circ(_a, getFinalFrequency());
}

double InspiralContainer::getTimeSpacing(){
  return _dt;
}
double InspiralContainer::getDuration(){
	return (getSize() - 1)*_dt;
}

int InspiralContainer::getSize(){
  return _alpha.size();
}

InspiralGenerator::InspiralGenerator(TrajectorySpline2D &traj, int num_threads): _traj(traj){
  if(num_threads > 0){
    omp_set_num_threads(num_threads);
  }
}

InspiralContainer InspiralGenerator::computeInspiral(double a, double massratio, double r0, double dt, double T, int num_threads){
	double chi, omega_i, alpha_i, t_i;
	computeInitialConditions(chi, omega_i, alpha_i, t_i, a, massratio, r0, T);
	int timesteps = computeTimeStepNumber(dt, T);
	InspiralContainer inspiral(timesteps);
	inspiral.setInspiralInitialConditions(a, massratio, r0, dt);
	computeInspiral(inspiral, chi, omega_i, alpha_i, t_i, massratio, dt, num_threads);
	return inspiral;
}

void InspiralGenerator::computeInitialConditions(double &chi, double &omega_i, double &alpha_i, double &t_i, double a, double massratio, double r0, double &T){
	chi = chi_of_spin(a);
	omega_i = kerr_geo_azimuthal_frequency_circ_time(a, r0);
	alpha_i = alpha_of_a_omega(a, omega_i, _traj.orbital_frequency_isco(chi));
	t_i = _traj.time(chi, alpha_i);
	if(T > -t_i/massratio){
		T = -t_i/massratio;
	}
}

void InspiralGenerator::computeInspiral(InspiralContainer &inspiral, double chi, double omega_i, double alpha_i, double t_i, double massratio, double dt, int num_threads){
	int steps = inspiral.getSize();
	double a = inspiral.getSpin();
	dt *= massratio; // need to rescale by massratio to get in terms of "slow time"

	double phase_i = _traj.phase_of_time(chi, t_i);
	inspiral.setTimeStep(0, alpha_i, 0.);	

	// first we downsample by fixing a and sampling in t
	// int downsample_steps = 200;
	// double t_f = t_i + dt*(time_steps-1);
	// double dt_downsample = (t_f - t_i)/double(downsample_steps - 1.);
	// Vector alpha_downsample(downsample_steps);
	// Vector phase_downsample(downsample_steps);
	// Vector t_downsample(downsample_steps);

	// #pragma omp parallel num_threads(num_threads)
	// {
	// 	#pragma omp for
	// 	for(int j = 0; j < downsample_steps - 1; j++){
	// 		t_downsample[j] = t_i + dt_downsample*j;
	// 		alpha_downsample[j] = _traj.orbital_alpha(chi, t_downsample[j]);
	// 		phase_downsample[j] = (_traj.phase(chi, alpha_downsample[j]) - phase_i)/massratio;
	// 	}
	// }

	// t_downsample[downsample_steps - 1] = t_f;
	// alpha_downsample[downsample_steps - 1] = _traj.orbital_alpha(chi, t_downsample[downsample_steps - 1]);
	// phase_downsample[downsample_steps - 1] = (_traj.phase(chi, alpha_downsample[downsample_steps - 1]) - phase_i)/massratio;
	// CubicSpline alpha_interp(t_downsample, alpha_downsample);
	// CubicSpline phase_interp(t_downsample, phase_downsample);

	// #pragma omp parallel num_threads(num_threads)
	// {
	// 	double alpha, phase, t_j;
	// 	#pragma omp for
	// 	for(int j = 0; j < time_steps; j++){
	// 		t_j = t_i + j*dt;
	// 		alpha = alpha_interp.evaluate(t_j);
	// 		phase = phase_interp.evaluate(t_j);
	// 		inspiral.setTimeStep(j, alpha, phase);
	// 	}
	// }
	#pragma omp parallel num_threads(num_threads)
	{
		double alpha, phase;
		#pragma omp for
		for(int j = 1; j < steps; j++){
			alpha = _traj.orbital_alpha(chi, t_i + dt*j);
			phase = (_traj.phase_of_time(chi, t_i + dt*j) - phase_i)/massratio;
			inspiral.setTimeStep(j, alpha, phase);
		}
	}
}

double InspiralGenerator::computeTimeToMerger(double a, double massratio, double r0){
	double chi = chi_of_spin(a);
	double omega_i = kerr_geo_azimuthal_frequency_circ_time(a, r0);
	double alpha_i = alpha_of_a_omega(a, omega_i, _traj.orbital_frequency_isco(chi));
	return -_traj.time(chi, alpha_i)/massratio;
}

int InspiralGenerator::computeTimeStepNumber(double a, double massratio, double r0, double dt, double T){
	double Tmerge = computeTimeToMerger(a, massratio, r0);
	T = (T > Tmerge) ? Tmerge : T;
	return computeTimeStepNumber(dt, T);
}

int InspiralGenerator::computeTimeStepNumber(double dt, double T){
	return T/dt + 1;
}

TrajectorySpline2D& InspiralGenerator::getTrajectorySpline(){
	return _traj;
}

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

	Data data = {
		.x = xVec,
		.y = yVec,
		.z = zVec
	};

	return data;
}

inline bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

TrajectoryData::TrajectoryData(const Vector &chi, const Vector &alpha, const Vector &t, const Vector &phi, const Vector &chiFlux, const Vector &alphaFlux, const Vector & flux, const Vector &beta, const Vector &omega, const Vector &alphaOfT, const Vector &phiOfT, const double &tMax): chi(chi), alpha(alpha), t(t), phi(phi), chiFlux(chiFlux), alphaFlux(alphaFlux), flux(flux), beta(beta), omega(omega), alphaOfT(alphaOfT), phiOfT(phiOfT), tMax(tMax) {}

TrajectoryData read_trajectory_data(std::string filename){
	int chiSample, alphaSample, chiSampleFlux, alphaSampleFlux;
	double chi, alpha, t, phi, flux, beta, omega, PhiT;
	
	if(!file_exists(filename)){
		return TrajectoryData();
	}

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
	lin >> chiSample >> alphaSample >> chiSampleFlux >> alphaSampleFlux;

	int n = chiSample*alphaSample;
	int nflux = chiSampleFlux*alphaSampleFlux;
	if(n > 500000){
		n = 1;
		std::cout << "ERROR: File "<< filename << " does not have appropriate number of samples \n";
	}
	Vector chiA(n), alphaA(n), tA(n), phiA(n), betaA(n), omegaA(n), alphaOfTA(n), phiOfTA(n);
	Vector chiAFlux(nflux), alphaAFlux(nflux), fluxA(nflux);
	int i = 0;
	int j = 0;

	double oISCO, freq, spin;
	for(std::string line; std::getline(inFile, line);){
		lin.clear();
		lin.str(line);
		if(lin >> chi >> alpha >> flux >> t >> phi >> beta >> omega >> PhiT){
			oISCO = kerr_isco_frequency(spin_of_chi(chi));
			chiA[i] = chi;
			alphaA[i] = alpha;
			spin = spin_of_chi(chi);
			freq = omega_of_a_alpha(spin, alpha);
			tA[i] = -t/normalize_time(freq, oISCO);
			phiA[i] = -phi/normalize_phase(freq, oISCO);
			if(flux > 0){
				if(j > nflux - 1){
					std::cout << "(TRAJ) Error: More flux samples than expected\n";
				}
				fluxA[j] = flux/normalize_energy_flux(freq);
				chiAFlux[j] = chi;
				alphaAFlux[j] = alpha;
				j++;			
			}
			betaA[i] = beta;
			omegaA[i] = omega;
			alphaOfTA[i] = alpha_of_a_omega(spin, omega, oISCO);
			phiOfTA[i] = log1p(-PhiT);
			i++;
		}
	}
	Vector alphaAReduce(alphaSample);
	for(int i = 0; i < alphaSample; i++){
		alphaAReduce[i] = alphaA[i];
	}
	Vector chiAReduce(chiSample);
	for(int i = 0; i < chiSample; i++){
		chiAReduce[i] = chiA[i*alphaSample];
	}
	Vector alphaAFluxReduce(alphaSampleFlux);
	for(int i = 0; i < alphaSampleFlux; i++){
		alphaAFluxReduce[i] = alphaAFlux[i];
	}
	Vector chiAFluxReduce(chiSampleFlux);
	for(int i = 0; i < chiSampleFlux; i++){
		chiAFluxReduce[i] = chiAFlux[i*alphaSample];
	}
	Vector betaAReduce(alphaSample);
	for(int i = 0; i < alphaSample; i++){
		betaAReduce[i] = betaA[i];
	}
	// double tMax = tA[n - 1];
	// double tMax = -t;
	double chiTemp = chiA[n - 1];
	double alphaTemp = alphaA[n - 1];
	double aTemp = spin_of_chi(chiTemp);
	double oISCOTemp = kerr_isco_frequency(spin_of_chi(chiTemp));
	double omegaTemp = omega_of_a_alpha(aTemp, alphaTemp);
	double tMax = tA[n - 1]*normalize_time(omegaTemp, oISCOTemp);

	TrajectoryData traj(chiAReduce, alphaAReduce, tA, phiA, chiAFluxReduce, alphaAFluxReduce, fluxA, betaAReduce, omegaA, alphaOfTA, phiOfTA, tMax);

	return traj;
}

CubicSpline read_time_spline(double){
  TrajectoryData traj = read_trajectory_data();
  return CubicSpline(traj.alpha, traj.t);
}

CubicSpline read_phase_spline(double){
  TrajectoryData traj = read_trajectory_data();
  return CubicSpline(traj.alpha, traj.phi);
}

CubicSpline read_flux_spline(double){
  TrajectoryData traj = read_trajectory_data();
  return CubicSpline(traj.alpha, traj.flux);
}

double slow_time_of_omega(CubicSpline &t, double a, double omega){
  return -t.evaluate(alpha_of_a_omega(a, omega));
}

double slow_phase_of_omega(CubicSpline &phi, double a, double omega){
  return -phi.evaluate(alpha_of_a_omega(a, omega));
}

double normalized_flux_of_omega(CubicSpline &Edot, double a, double omega){
  return Edot.evaluate(alpha_of_a_omega(a, omega))*newtonian_energy_flux(omega);
}

double normalized_omega_time_derivative_of_omega(CubicSpline &Edot, double a, double omega){
  return -(1./kerr_geo_denergy_domega_circ(a, omega) + omega/kerr_geo_dmomentum_domega_circ(a, omega))*Edot.evaluate(alpha_of_a_omega(a, omega))*newtonian_energy_flux(omega);
}

// TrajectorySpline2D class

TrajectorySpline2D::TrajectorySpline2D(std::string filename): TrajectorySpline2D(read_trajectory_data(filename)) {}
TrajectorySpline2D::TrajectorySpline2D(TrajectoryData traj):
TrajectorySpline2D(traj.chi, traj.alpha, traj.chiFlux, traj.alphaFlux, traj.beta, traj.t, traj.phi, traj.flux, traj.omega, traj.alphaOfT, traj.phiOfT, traj.tMax) {}
TrajectorySpline2D::TrajectorySpline2D(const Vector & chi, const Vector & alpha, const Vector & chiFlux, const Vector & alphaFlux, const Vector & beta, const Vector & t, const Vector & phi, const Vector & flux, const Vector & omega, const Vector & alphaOfT, const Vector & phaseOfT, const double & tMax):
  	_time_spline(chi, alpha, t), _phase_spline(chi, alpha, phi), _flux_spline(chiFlux, alphaFlux, flux), _alpha_spline(chi, beta, alphaOfT), _frequency_spline(chi, beta, omega), _phase_time_spline(chi, beta, phaseOfT), _time_norm_parameter(gamma_of_time(-tMax)) {}
TrajectorySpline2D::~TrajectorySpline2D(){}

// Frequency domain

double TrajectorySpline2D::time(double chi, double alpha){
  	// return -expm1(pow(_time_spline.evaluate(chi, alpha), 2));
	double oISCO = abs(kerr_isco_frequency(spin_of_chi(chi)));
	double omega = omega_of_a_alpha(chi, alpha, oISCO);
	return -_time_spline.evaluate(chi, alpha)*normalize_time(omega, oISCO);
}

double TrajectorySpline2D::phase(double chi, double alpha){
  	// return -expm1(pow(_phase_spline.evaluate(chi, alpha), 2));
	double oISCO = kerr_isco_frequency(spin_of_chi(chi));
	double omega = omega_of_a_alpha(chi, alpha, oISCO);
	return -_phase_spline.evaluate(chi, alpha)*normalize_phase(omega, oISCO);
}

double TrajectorySpline2D::flux(double chi, double alpha){
	double omega = omega_of_a_alpha(spin_of_chi(chi), alpha);
  	return _flux_spline.evaluate(chi, alpha)*normalize_energy_flux(omega);
}

double TrajectorySpline2D::flux_norm(double chi, double alpha){
	return _flux_spline.evaluate(chi, alpha);
}

double TrajectorySpline2D::orbital_frequency_time_derivative_from_flux(double chi, double alpha){
	double a = spin_of_chi(chi);
	double omega = omega_of_a_alpha(a, alpha);
  	return -(1./kerr_geo_denergy_domega_circ(a, omega))*_flux_spline.evaluate(chi, alpha)*normalize_energy_flux(omega);
}

double TrajectorySpline2D::time_of_a_omega(double a, double omega){
  	return time(chi_of_spin(a), alpha_of_a_omega(a, omega));
}

double TrajectorySpline2D::time_of_a_omega_derivative(double a, double omega){
	double oISCO = abs(kerr_isco_frequency(a));
	double chi = chi_of_spin(a);
	double alpha = alpha_of_a_omega(a, omega, oISCO);
	// double f = _time_spline.evaluate(chi, alpha);
  	// return -2.*f*exp(pow(f, 2))*_time_spline.derivative_y(chi, alpha)*dalpha_domega_of_a_omega(a, omega, oISCO);
	return _time_spline.evaluate(chi, alpha)*normalize_time_domega(omega) + dalpha_domega_of_a_omega(a, omega, oISCO)*_time_spline.derivative_y(chi, alpha)*normalize_time(omega, oISCO);
}

// double TrajectorySpline2D::time_of_a_alpha_omega_derivative(double a, double alpha){
// 	double oISCO = abs(kerr_isco_frequency(a));
// 	double chi = chi_of_spin(a);
// 	// double f = _time_spline.evaluate(chi, alpha);
//   	// return -2.*f*exp(pow(f, 2))*_time_spline.derivative_y(chi, alpha)*dalpha_domega_of_alpha(alpha, oISCO);
// 	double omega = omega_of_a_alpha(a, alpha, oISCO);
// 	return _time_spline.evaluate(chi, alpha)*normalize_time_domega(omega) + dalpha_domega_of_a_omega(a, omega, oISCO)*_time_spline.derivative_y(chi, alpha)*normalize_time(omega, oISCO);
// }

double TrajectorySpline2D::time_of_a_alpha_omega_derivative(double a, double alpha){
	double oISCO = abs(kerr_isco_frequency(a));
	double chi = chi_of_spin(a);
	// double f = _time_spline.evaluate(chi, alpha);
  	// return -2.*f*exp(pow(f, 2))*_time_spline.derivative_y(chi, alpha)*dalpha_domega_of_alpha(alpha, oISCO);
	double omega = omega_of_a_alpha(a, alpha, oISCO);
	return -(kerr_geo_denergy_domega_circ(a, omega))/_flux_spline.evaluate(chi, alpha)/normalize_energy_flux(omega);
}

double TrajectorySpline2D::phase_of_a_omega(double a, double omega){
  	return phase(chi_of_spin(a), alpha_of_a_omega(a, omega));
}

double TrajectorySpline2D::phase_of_a_omega_derivative(double a, double omega){
	double oISCO = abs(kerr_isco_frequency(a));
	double chi = chi_of_spin(a);
	double alpha = alpha_of_a_omega(a, omega, oISCO);
	// double f = _phase_spline.evaluate(chi, alpha);
  	// return -2.*f*exp(pow(f, 2))*_phase_spline.derivative_y(chi, alpha)*dalpha_domega_of_a_omega(a, omega, oISCO);
	return _phase_spline.evaluate(chi, alpha)*normalize_phase_domega(omega) + dalpha_domega_of_a_omega(a, omega, oISCO)*_phase_spline.derivative_y(chi, alpha)*normalize_phase(omega, oISCO);
}

double TrajectorySpline2D::flux_of_a_omega(double a, double omega){
  	return _flux_spline.evaluate(chi_of_spin(a), alpha_of_a_omega(a, omega))*normalize_energy_flux(omega);
}

double TrajectorySpline2D::orbital_frequency_time_derivative_from_flux_of_a_omega(double a, double omega){
  	return -(1./kerr_geo_denergy_domega_circ(a, omega))*_flux_spline.evaluate(chi_of_spin(a), alpha_of_a_omega(a, omega))*normalize_energy_flux(omega);
}

// Time domain

double TrajectorySpline2D::orbital_alpha(double chi, double t){
  	return _alpha_spline.evaluate(chi, beta_of_time(t, _time_norm_parameter));
}

double TrajectorySpline2D::orbital_alpha_derivative(double chi, double t){
  	return _alpha_spline.derivative_y(chi, beta_of_time(t, _time_norm_parameter))*dbeta_dtime(t, _time_norm_parameter);
}

double TrajectorySpline2D::phase_of_time(double chi, double t){
  	return phase_of_phase_norm_2(_phase_time_spline.evaluate(chi, beta_of_time(t, _time_norm_parameter)));
}

double TrajectorySpline2D::phase_of_time_derivative(double chi, double t){
	double phase2 = _phase_time_spline.derivative_y(chi, beta_of_time(t, _time_norm_parameter));
  	return _phase_time_spline.derivative_y(chi, beta_of_time(t, _time_norm_parameter))*dbeta_dtime(t, _time_norm_parameter)*dphase_dphase_norm_2(phase2);
}

double TrajectorySpline2D::phase_of_a_time(double a, double t){
	double chi = chi_of_spin(a);
  	return phase_of_phase_norm_2(_phase_time_spline.evaluate(chi, beta_of_time(t, _time_norm_parameter)));
}

double TrajectorySpline2D::phase_of_a_time_derivative(double a, double t){
	double chi = chi_of_spin(a);
	double phase2 = _phase_time_spline.derivative_y(chi, beta_of_time(t, _time_norm_parameter));
  	return _phase_time_spline.derivative_y(chi, beta_of_time(t, _time_norm_parameter))*dbeta_dtime(t, _time_norm_parameter)*dphase_dphase_norm_2(phase2);
}

double TrajectorySpline2D::orbital_frequency(double a, double t){
  	return _frequency_spline.evaluate(chi_of_spin(a), beta_of_time(t, _time_norm_parameter));
}

double TrajectorySpline2D::orbital_frequency_derivative(double a, double t){
  	return _frequency_spline.derivative_y(chi_of_spin(a), beta_of_time(t, _time_norm_parameter))*dbeta_dtime(t, _time_norm_parameter);
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
	// return 1. - exp(pow(_max_time_spline.evaluate(chi_of_spin(a)), 2));
	return _time_spline.evaluate(chi_of_spin(a), ALPHA_MAX);
}

void TrajectorySpline2D::flux_of_a_omega(double flux[], const double a[], const double omega[], int n, int num_threads){
	if(num_threads > 1){
        omp_set_num_threads(num_threads);
    }
	#pragma omp parallel for
		for(int j = 0; j < n; j++){
			flux[j] = _flux_spline.evaluate(chi_of_spin(a[j]), alpha_of_a_omega(a[j], omega[j]))*normalize_energy_flux(omega[j]);
		}
}