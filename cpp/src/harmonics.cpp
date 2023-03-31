#include "harmonics.hpp"

HarmonicModeData read_harmonic_mode_data(int L, int m, std::string filepath_base){
	std::string filepath = filepath_base + "_" + std::to_string(L) + "_" + std::to_string(m) + ".txt";
	double chi, alpha, Alm, Philm;
	Vector chiVec, alphaVec, AlmVec, PhilmVec;
	std::istringstream lin;
	std::ifstream inFile(filepath);
	std::string line;
	std::getline(inFile, line);
	lin.clear();
	lin.str(line);
	int chiSample, alphaSample;
	while(line.front() == '#' || line.empty() || isalpha(line.front())){
		std::getline(inFile, line);
		lin.clear();
		lin.str(line);
	}
	lin >> chiSample >> alphaSample;
	int n = chiSample*alphaSample;
	EigenArray chiA(n), alphaA(n), AlmA(n), PhilmA(n);
	int i = 0;

	for(std::string line; std::getline(inFile, line);){
		lin.clear();
		lin.str(line);
		if(lin >> chi >> alpha >> Alm >> Philm){
			chiA[i] = chi;
			alphaA[i] = alpha;
			AlmA[i] = Alm;
			PhilmA[i] = Philm;
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

	HarmonicModeData mode = {
    	.chi = chiAReduce,
		.alpha = alphaAReduce,
		.A = AlmA,
		.Phi = PhilmA
	};

	return mode;
}

HarmonicSpline::HarmonicSpline(double chi, const EigenArray & alpha, const EigenArray & A, const EigenArray & Phi): _spin(spin_of_chi(chi)), _amplitude_spline(alpha, A), _phase_spline(alpha, Phi) {}
HarmonicSpline::HarmonicSpline(double spin, EigenCubicInterpolator amplitude_spline, EigenCubicInterpolator phase_spline): _spin(spin), _amplitude_spline(amplitude_spline), _phase_spline(phase_spline) {}
HarmonicSpline::~HarmonicSpline() {}

double HarmonicSpline::amplitude(double alpha){
  return _amplitude_spline.evaluate(alpha);
}

double HarmonicSpline::phase(double alpha){
  return _phase_spline.evaluate(alpha);
}

double HarmonicSpline::amplitude_of_omega(double omega){
  return _amplitude_spline.evaluate(alpha_of_a_omega(_spin, omega));
}

double HarmonicSpline::phase_of_omega(double omega){
  return _phase_spline.evaluate(alpha_of_a_omega(_spin, omega));
}

double HarmonicSpline::phase_of_omega_derivative(double omega){
  return _phase_spline.derivative(alpha_of_a_omega(_spin, omega))*dalpha_domega_of_a_omega(_spin, omega, abs(kerr_isco_frequency(_spin)));
}

HarmonicSpline2D::HarmonicSpline2D(int L, int m, std::string filepath_base): HarmonicSpline2D(read_harmonic_mode_data(L, m, filepath_base)) {}
HarmonicSpline2D::HarmonicSpline2D(HarmonicModeData mode): _amplitude_spline(mode.chi, mode.alpha, mode.A), _phase_spline(mode.chi, mode.alpha, mode.Phi) { }
HarmonicSpline2D::HarmonicSpline2D(const EigenArray & chi, const EigenArray & alpha, const EigenArray & Amp, const EigenArray & Phi): _amplitude_spline(chi, alpha, Amp), _phase_spline(chi, alpha, Phi) { }
HarmonicSpline2D::~HarmonicSpline2D() {}

double HarmonicSpline2D::amplitude(double chi, double alpha){
  return _amplitude_spline.evaluate(chi, alpha);
}

double HarmonicSpline2D::phase(double chi, double alpha){
  return _phase_spline.evaluate(chi, alpha);
}

double HarmonicSpline2D::amplitude_of_a_omega(double a, double omega){
  return _amplitude_spline.evaluate(chi_of_spin(a), alpha_of_a_omega(a, omega));
}

double HarmonicSpline2D::phase_of_a_omega(double a, double omega){
  return _phase_spline.evaluate(chi_of_spin(a), alpha_of_a_omega(a, omega));
}

double HarmonicSpline2D::phase_of_a_omega_derivative(double a, double omega){
  return _phase_spline.derivative_y(chi_of_spin(a), alpha_of_a_omega(a, omega))*dalpha_domega_of_a_omega(a, omega, abs(kerr_isco_frequency(a)));
}

// EigenBicubicInterpolator HarmonicSpline2D::getAmplitudeSpline(){
// 	return _amplitude_spline;
// }
// EigenBicubicInterpolator HarmonicSpline2D::getPhaseSpline(){
// 	return _phase_spline;
// }
// EigenCubicInterpolator HarmonicSpline2D::getReducedAmplitudeSpline(double a){
// 	return _amplitude_spline.reduce_x(chi_of_spin(a));
// }
// EigenCubicInterpolator HarmonicSpline2D::getReducedPhaseSpline(double a){
// 	return _phase_spline.reduce_x(chi_of_spin(a));
// }

// HarmonicSpline HarmonicSpline2D::constant_spin_spline(double a){
// 	double chi = chi_of_spin(a);
// 	return HarmonicSpline(a, _amplitude_spline.reduce_x(chi), _phase_spline.reduce_x(chi));
// }

///////////////////////////////////////////////////////
//////////        HarmonicAmplitudes     //////////////
///////////////////////////////////////////////////////

// A Harmonic class that holds several different modes
HarmonicAmplitudes::HarmonicAmplitudes(std::vector<int> &lmodes, std::vector<int> &mmodes, std::string filepath_base): 
	HarmonicAmplitudes(lmodes.data(), mmodes.data(), lmodes.size(), filepath_base) {}

HarmonicAmplitudes::HarmonicAmplitudes(int lmodes[], int mmodes[], int modeNum, std::string filepath_base): _modeNum(modeNum), _harmonics(modeNum) {
	for(int i = 0; i < _modeNum; i++){
		_harmonics[i] = new HarmonicSpline2D(lmodes[i], mmodes[i], filepath_base);
		_position_map[std::pair<int, int>(lmodes[i], mmodes[i])] = i;
	}
}
HarmonicAmplitudes::~HarmonicAmplitudes() {
	for(int i = 0; i < _modeNum; i++){
		delete _harmonics[i];
	}
}

int HarmonicAmplitudes::key_check(std::pair<int, int> key){
	if(_position_map.count(key) > 0){
		return 1;
	}else{
		return 0;
	}
}

double HarmonicAmplitudes::amplitude(int l, int m, double chi, double alpha){
	std::pair<int, int> key(l, m);
	return _harmonics[_position_map[key]]->amplitude(chi, alpha);
}

double HarmonicAmplitudes::phase(int l, int m, double chi, double alpha){
	std::pair<int, int> key(l, m);
	return _harmonics[_position_map[key]]->phase(chi, alpha);
}

double HarmonicAmplitudes::amplitude_of_a_omega(int l, int m, double a, double omega){
	std::pair<int, int> key(l, m);
	return _harmonics[_position_map[key]]->amplitude_of_a_omega(a, omega);
}

double HarmonicAmplitudes::phase_of_a_omega(int l, int m, double a, double omega){
	std::pair<int, int> key(l, m);
	return _harmonics[_position_map[key]]->phase_of_a_omega(a, omega);
}

double HarmonicAmplitudes::phase_of_a_omega_derivative(int l, int m, double a, double omega){
	std::pair<int, int> key(l, m);
	return _harmonics[_position_map[key]]->phase_of_a_omega_derivative(a, omega);
}

HarmonicSpline2D* HarmonicAmplitudes::getPointer(int l, int m){
	std::pair<int, int> key(l, m);
	return _harmonics[_position_map[key]];
}

// double HarmonicAmplitudes::amplitude(int l, int m, double chi, double alpha){
// 	std::pair<int, int> key(l, m);
// 	if(key_check(key)){
// 		return _harmonics[_position_map[key]]->amplitude(chi, alpha);
// 	}else{
// 		std::cout << "Invalid key\n";
// 		return 0.;
// 	}
// }

// double HarmonicAmplitudes::phase(int l, int m, double chi, double alpha){
// 	std::pair<int, int> key(l, m);
// 	if(key_check(key)){
// 		return _harmonics[_position_map[key]]->phase(chi, alpha);
// 	}else{
// 		std::cout << "Invalid key\n";
// 		return 0.;
// 	}
// }

// double HarmonicAmplitudes::amplitude_of_a_omega(int l, int m, double a, double omega){
// 	std::pair<int, int> key(l, m);
// 	if(key_check(key)){
// 		return _harmonics[_position_map[key]]->amplitude_of_a_omega(a, omega);
// 	}else{
// 		std::cout << "Invalid key\n";
// 		return 0.;
// 	}
// }

// double HarmonicAmplitudes::phase_of_a_omega(int l, int m, double a, double omega){
// 	std::pair<int, int> key(l, m);
// 	if(key_check(key)){
// 		return _harmonics[_position_map[key]]->phase_of_a_omega(a, omega);
// 	}else{
// 		std::cout << "Invalid key\n";
// 		return 0.;
// 	}
// }

// double HarmonicAmplitudes::phase_of_a_omega_derivative(int l, int m, double a, double omega){
// 	std::pair<int, int> key(l, m);
// 	if(key_check(key)){
// 		return _harmonics[_position_map[key]]->phase_of_a_omega_derivative(a, omega);
// 	}else{
// 		std::cout << "Invalid key\n";
// 		return 0.;
// 	}
// }

HarmonicSelector::HarmonicSelector(HarmonicAmplitudes &harm, HarmonicOptions opts): _harm(harm), _opts(opts) {}

double HarmonicSelector::modePower(int l, int m, InspiralContainer &inspiral){
	return modePower(l, m, inspiral, _opts);
}
int HarmonicSelector::gradeMode(int l, int m, InspiralContainer &inspiral, double power22){
	return gradeMode(l, m, inspiral, power22, _opts);
}
int HarmonicSelector::gradeMode(int l, int m, InspiralContainer &inspiral, double power22, double &plusYlm, double &crossYlm, double theta){
	return gradeMode(l, m, inspiral, power22, plusYlm, crossYlm, theta, _opts);
}
HarmonicModeContainer HarmonicSelector::selectModes(InspiralContainer &inspiral, double theta){
	return selectModes(inspiral, theta, _opts);
}

double HarmonicSelector::modePower(int l, int m, InspiralContainer &inspiral, HarmonicOptions opts){
	double power = 0.;
	int timeSteps = inspiral.getSize();
	int stepSize = 1;
	int max_samples = opts.max_samples;
	double chi = chi_of_spin(inspiral.getSpin());
	if(timeSteps > max_samples){
		stepSize = int(timeSteps/max_samples);
	}else{
		max_samples = timeSteps;
	}
	double amp2;
	for(int i = 0; i < max_samples; i++){
		amp2 = pow(_harm.amplitude(l, m, chi, inspiral.getAlpha(timeSteps - 1 - i*stepSize)), 2);
		power += _harm.amplitude(l, m, chi, inspiral.getAlpha(timeSteps - 1 - i*stepSize));
	}

	return power;
}

int HarmonicSelector::gradeMode(int l, int m, InspiralContainer &inspiral, double power22, HarmonicOptions opts){
	double powerLM = modePower(l, m, inspiral, opts);
	if(powerLM/power22 > opts.epsilon){
		return 1;
	}else{
		return 0;
	}
}

int HarmonicSelector::gradeMode(int l, int m, InspiralContainer &inspiral, double power22, double &plusYlm, double &crossYlm, double theta, HarmonicOptions opts){
	double powerLM = modePower(l, m, inspiral, opts);
	Yslm_plus_cross_polarization(plusYlm, crossYlm, l, m, theta);
	if(powerLM*(pow(plusYlm, 2) + pow(crossYlm, 2))/power22 > opts.epsilon){
		return 1;
	}else{
		return 0;
	}
}

HarmonicModeContainer HarmonicSelector::selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts){
	HarmonicModeContainer harmonics;
	double plusY, crossY, power22;
	int l, m;
	
	l = 2;
	m = 2;
	power22 = modePower(l, m, inspiral, opts);
	Yslm_plus_cross_polarization(plusY, crossY, l, m, theta);
	power22 *= pow(plusY, 2) + pow(crossY, 2);

	harmonics.lmodes.push_back(l);
	harmonics.mmodes.push_back(m);
	harmonics.plusY.push_back(plusY);
	harmonics.crossY.push_back(crossY);

	m = 1;
	if(gradeMode(l, m, inspiral, power22, plusY, crossY, theta, opts)){
		harmonics.lmodes.push_back(l);
		harmonics.mmodes.push_back(m);
		harmonics.plusY.push_back(plusY);
		harmonics.crossY.push_back(crossY);
	}

	int lmax = 15;;
	for(l = 3; l < lmax; l++){
		for(m = l; m > 0; m--){
			if(gradeMode(l, m, inspiral, power22, plusY, crossY, theta, opts)){
				harmonics.lmodes.push_back(l);
				harmonics.mmodes.push_back(m);
				harmonics.plusY.push_back(plusY);
				harmonics.crossY.push_back(crossY);
			}else{
				m = 0;
			}
		}
	}

	return harmonics;
}

HarmonicOptions HarmonicSelector::getHarmonicOptions(){
	return _opts;
}


double Yslm_plus_polarization(int l, int m, double theta){
	double sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
    double sYlmMinus = spin_weighted_spherical_harmonic(2, l, m, theta);
	return sYlm + pow(-1, l+m)*sYlmMinus;
}
double Yslm_cross_polarization(int l, int m, double theta){
	double sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
    double sYlmMinus = spin_weighted_spherical_harmonic(2, l, m, theta);
	return sYlm - pow(-1, l+m)*sYlmMinus;
}
void Yslm_plus_cross_polarization(double &plusY, double &crossY, int l, int m, double theta){
	double sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
    double sYlmMinus = spin_weighted_spherical_harmonic(2, l, m, theta);
	plusY = sYlm + pow(-1, l+m)*sYlmMinus;
	crossY = sYlm - pow(-1, l+m)*sYlmMinus;
}