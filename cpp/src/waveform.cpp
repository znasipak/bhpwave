#include "waveform.hpp"

#define REL_MODE_EPS 5.e-4
#define OMEGA_DIFF_EPS 1.e-5
#define ALPHA_MAX 1.
#define ALPHA_MIN 0.

WaveformContainer::WaveformContainer(int timeSteps): _size(timeSteps), _owner_flag(1) {
	_plus = new double[timeSteps];
	_cross = new double[timeSteps];
  for(int i = 0; i < timeSteps; i++){
    _plus[i] = 0.;
    _cross[i] = 0.;
  }
}

WaveformContainer::WaveformContainer(double *plus_ptr, double *cross_ptr, int timeSteps): _plus(plus_ptr), _cross(cross_ptr), _size(timeSteps), _owner_flag(0) {}

WaveformContainer::~WaveformContainer(){
	if(_owner_flag){
		delete _plus;
		delete _cross;
	}
}

void WaveformContainer::setTimeStep(int i, double plus, double cross){
  _plus[i] = plus;
  _cross[i] = cross;
}

void WaveformContainer::addTimeStep(int i, double plus, double cross){
  #pragma omp atomic
  _plus[i] += plus;

  #pragma omp atomic
  _cross[i] += cross;
}

void WaveformContainer::multiplyTimeStep(int i, double plus, double cross){
  _plus[i] *= plus;
  _cross[i] *= cross;
}

double* WaveformContainer::getPlusPointer(){
	return _plus;
}
double* WaveformContainer::getCrossPointer(){
	return _cross;
}

// Vector& WaveformContainer::getPlusPolarizationNonConstRef(){
//     return _plus;
// }
// Vector& WaveformContainer::getCrossPolarizationNonConstRef(){
//     return _cross;
// }

double WaveformContainer::getPlus(int i){
  return _plus[i];
}

double WaveformContainer::getCross(int i){
  return _cross[i];
}

int WaveformContainer::getSize(){
  return _size;
}

WaveformHarmonicGenerator::WaveformHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts): _Alm(Alm), _mode_selector(Alm, hOpts), _opts(wOpts) {}

WaveformContainer WaveformHarmonicGenerator::computeWaveformHarmonic(int l, int m, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts){
	WaveformContainer h(inspiral.getSize());
	computeWaveformHarmonic(h, l, m, inspiral, theta, phi, opts);
	return h;
}
WaveformContainer WaveformHarmonicGenerator::computeWaveformHarmonics(int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts){
	WaveformContainer h(inspiral.getSize());
	computeWaveformHarmonics(h, l, m, modeNum, inspiral, theta, phi, opts);
	return h;
}

void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi){
	computeWaveformHarmonics(h, inspiral, theta, phi, _opts);
}
void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, HarmonicOptions hOpts){
	computeWaveformHarmonics(h, inspiral, theta, phi, hOpts, _opts);
}
void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions wOpts){
	computeWaveformHarmonics(h, inspiral, theta, phi, getHarmonicOptions());
}
void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, InspiralContainer &inspiral, double theta, double phi, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	HarmonicModeContainer modes = _mode_selector.selectModes(inspiral, theta, hOpts);
	computeWaveformHarmonics(h, modes.lmodes.data(), modes.mmodes.data(), modes.plusY.data(), modes.crossY.data(), modes.lmodes.size(), inspiral, theta, phi, wOpts);
}


void WaveformHarmonicGenerator::computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi){
	computeWaveformHarmonic(h, l, m, inspiral, theta, phi, _opts);
}
void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi){
	computeWaveformHarmonics(h, l, m, modeNum, inspiral, theta, phi, _opts);
}
void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi){
	computeWaveformHarmonics(h, l, m, plusY, crossY, modeNum, inspiral, theta, phi, _opts);
}

void WaveformHarmonicGenerator::computeWaveformHarmonic(WaveformContainer &h, int l, int m, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts){    
    double sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
    double sYlmMinus = spin_weighted_spherical_harmonic(2, l, m, theta);
    double mphi_mod_2pi = fmod(m*phi, 2.*M_PI);
    double chi = chi_of_spin(inspiral.getSpin());
    double plusY = (sYlm + pow(-1, l + m)*sYlmMinus);
    double crossY = (sYlm - pow(-1, l + m)*sYlmMinus);

    int imax = inspiral.getSize();

    #pragma omp parallel num_threads(opts.num_threads)
    {
		double amp, modePhase, Phi, hplus, hcross;
        #pragma omp for
        for(int i = 0; i < imax; i++){
            amp = _Alm.amplitude(l, m, chi, inspiral.getAlpha(i));
            modePhase = _Alm.phase(l, m, chi, inspiral.getAlpha(i));
            Phi = modePhase - fmod(m*inspiral.getPhase(i), 2.*M_PI) + mphi_mod_2pi;
            hplus = amp*plusY*std::cos(Phi);
			      hcross = -amp*crossY*std::sin(Phi);
            h.addTimeStep(i, hplus, hcross);
        }
    }
}

void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts){
    double plusY[modeNum], crossY[modeNum], mphi_mod_2pi[modeNum];
    HarmonicSpline2D* Alms[modeNum];
    double twopi = 2.*M_PI;

    for(int i = 0; i < modeNum; i++){
      double sYlm = spin_weighted_spherical_harmonic(-2, l[i], m[i], theta);
      double sYlmMinus = pow(-1, l[i] + m[i])*spin_weighted_spherical_harmonic(2, l[i], m[i], theta);
      plusY[i] = (sYlm + sYlmMinus);
      crossY[i] = (sYlm - sYlmMinus);
      mphi_mod_2pi[i] = fmod(m[i]*phi, twopi);
      Alms[i] = _Alm.getPointer(l[i], m[i]);
    }
    double chi = chi_of_spin(inspiral.getSpin());

    const int imax = inspiral.getSize();
    Vector hplus(modeNum*imax);
    Vector hcross(modeNum*imax);

    #pragma omp parallel num_threads(opts.num_threads)
    {
      int i, j;
      double amp, modePhase, Phi;
      // first we calculate all of the mode data
      #pragma omp for collapse(2) schedule(static)
      for(j = 0; j < modeNum; j++){
        for(i = 0; i < imax; i++){
          amp = Alms[j]->amplitude(chi, inspiral.getAlpha(i));
          modePhase = Alms[j]->phase(chi, inspiral.getAlpha(i));
          Phi = modePhase - fmod(m[j]*inspiral.getPhase(i), twopi) + mphi_mod_2pi[j];
          hplus[j + i*modeNum] = amp*plusY[j]*std::cos(Phi);
          hcross[j + i*modeNum] = -amp*crossY[j]*std::sin(Phi);
        }
      }

      // we then sum over the data and have it stored in the waveform container
      #pragma omp for schedule(static)
      for(int i = 0; i < imax; i++){
        for(int j = 0; j < modeNum; j++){
          h.addTimeStep(i, hplus[j + i*modeNum], hcross[j + i*modeNum]);
        }
      }
    }
}

void WaveformHarmonicGenerator::computeWaveformHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, double theta, double phi, WaveformHarmonicOptions opts){
    double mphi_mod_2pi[modeNum];
    HarmonicSpline2D* Alms[modeNum];
    double twopi = 2.*M_PI;

    // first compute mode-dependent but not time-step dependent information and store
    for(int i = 0; i < modeNum; i++){
      mphi_mod_2pi[i] = fmod(m[i]*phi, twopi);
      Alms[i] = _Alm.getPointer(l[i], m[i]);
    }
    double chi = chi_of_spin(inspiral.getSpin());

    int imax = h.getSize();
    // for some reason creating arrays of this size was causing issues with OpenMP
    // minimal internet research suggests that OpenMP allocates arrays on the stack
    // and large arrays can lead to crashes. std::vector seems to work better here
    Vector hplus(modeNum*imax);
    Vector hcross(modeNum*imax);

    #pragma omp parallel num_threads(opts.num_threads)
    {
      int i, j;
      double amp, modePhase, Phi;
      // first we calculate all of the mode data
      #pragma omp for collapse(2) schedule(static)
      for(j = 0; j < modeNum; j++){
        for(i = 0; i < imax; i++){
          amp = Alms[j]->amplitude(chi, inspiral.getAlpha(i));
          modePhase = Alms[j]->phase(chi, inspiral.getAlpha(i));
          Phi = modePhase - fmod(m[j]*inspiral.getPhase(i), twopi) + mphi_mod_2pi[j];
          hplus[j + i*modeNum] = amp*plusY[j]*std::cos(Phi);
          hcross[j + i*modeNum] = -amp*crossY[j]*std::sin(Phi);
        }
      }

      // we then sum over the data and have it stored in the waveform container
      #pragma omp for schedule(static)
      for(int i = 0; i < imax; i++){
        for(int j = 0; j < modeNum; j++){
          h.addTimeStep(i, hplus[j + i*modeNum], hcross[j + i*modeNum]);
        }
      }
    }
}

HarmonicSelector& WaveformHarmonicGenerator::getModeSelector(){
	return _mode_selector;
}

HarmonicModeContainer WaveformHarmonicGenerator::selectModes(InspiralContainer &inspiral, double theta){
	return _mode_selector.selectModes(inspiral, theta);
}

HarmonicModeContainer WaveformHarmonicGenerator::selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts){
	return _mode_selector.selectModes(inspiral, theta, opts);
}

WaveformHarmonicOptions WaveformHarmonicGenerator::getWaveformHarmonicOptions(){
	return _opts;
}

HarmonicOptions WaveformHarmonicGenerator::getHarmonicOptions(){
	return _mode_selector.getHarmonicOptions();
}

// Waveform Generator

WaveformGenerator::WaveformGenerator(TrajectorySpline2D &traj, HarmonicAmplitudes &harm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts): WaveformHarmonicGenerator(harm, hOpts, wOpts), _inspiralGen(traj) {}

double WaveformGenerator::convertTime(double t, double M){
	return t/solar_mass_to_seconds(M);
}

int WaveformGenerator::computeTimeStepNumber(double dt, double T){
	return years_to_seconds(T)/dt + 1;
}

int WaveformGenerator::computeTimeStepNumber(double M, double mu, double a, double r0, double dt, double T){
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	return _inspiralGen.computeTimeStepNumber(a, mu/M, r0, dt, T);
}

void WaveformGenerator::computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T){
	computeWaveform(h, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, getHarmonicOptions(), _opts);
}

void WaveformGenerator::computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts){
	computeWaveform(h, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, hOpts, _opts);
}

void WaveformGenerator::computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, WaveformHarmonicOptions wOpts){
	computeWaveform(h, M, mu, a, r0, dist, qS, phiS, qK, phiK, Phi_phi0, dt, T, getHarmonicOptions(), wOpts);
}

void WaveformGenerator::computeWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_strain_amplitude(mu, dist);

	// omp_set_num_threads(16);
	StopWatch watch;
	// watch.start();
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, wOpts.num_threads);
	// watch.stop();
	// watch.print();
	// watch.reset();
	computeWaveformHarmonics(h, inspiral, theta, phi - Phi_phi0, hOpts, wOpts);
	
	double rescaleRe, rescaleIm;
	rescaleRe = std::real(wOpts.rescale);
	rescaleIm = std::imag(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes

	// watch.start();
	int imax = h.getSize();
	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		double hplus, hcross;
		#pragma omp for
		for(i = 0; i < imax; i++){
		// total_td_check = omp_get_num_threads();
			hplus = h.getPlus(i);
			hcross = h.getCross(i);
			h.setTimeStep(i, rescaleRe*hplus + rescaleIm*hcross, rescaleRe*hcross - rescaleIm*hplus);
		}
	}
	// watch.stop();
	// watch.print();
	// watch.reset();
}

void WaveformGenerator::computeWaveform(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
  	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_strain_amplitude(mu, dist);

	// omp_set_num_threads(16);
	// StopWatch watch;
	// watch.start();
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, wOpts.num_threads);
	// watch.stop();
	// watch.print();
	// watch.reset();
	computeWaveformHarmonics(h, l, m, modeNum, inspiral, theta, phi - Phi_phi0, wOpts);
	
	double rescaleRe, rescaleIm;
	rescaleRe = std::real(wOpts.rescale);
	rescaleIm = std::imag(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes

	// watch.start();
	int imax = h.getSize();
	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		double hplus, hcross;
		#pragma omp for
		for(i = 0; i < imax; i++){
		// total_td_check = omp_get_num_threads();
			hplus = h.getPlus(i);
			hcross = h.getCross(i);
			h.setTimeStep(i, rescaleRe*hplus + rescaleIm*hcross, rescaleRe*hcross - rescaleIm*hplus);
		}
	}
	// watch.stop();
	// watch.print();
	// watch.reset();
}

HarmonicModeContainer WaveformGenerator::selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T){
	return selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, getHarmonicOptions());
}
HarmonicModeContainer WaveformGenerator::selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions opts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions wOpts = getWaveformHarmonicOptions();

	dt = T/(opts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, wOpts.num_threads);
	return WaveformHarmonicGenerator::selectModes(inspiral, theta, opts);
}

void WaveformGenerator::computeWaveformSourceFrame(WaveformContainer &h, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T){
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions opts = getWaveformHarmonicOptions();

	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, opts.num_threads);
	computeWaveformHarmonics(h, inspiral, theta, phi - Phi_phi0, opts);
}

void sourceAngles(double &theta, double &phi, double qS, double phiS, double qK, double phiK){
    // Calculate the sky location :math:`(theta, phi)` of the observor in the
    // source frame using the sky location and orientation of the source in
    // the SSB frame

    // args:
    //     qS (double): polar angle of the source's sky location
    //     phiS (double): azimuthal angle of the source's sky location
    //     qK (double): polar angle of the Kerr spin vector
    //     phiK (double): azimuthal angle of the Kerr spin vector

    phi = -0.5*M_PI;
    theta = acos(-(sin(qS)*sin(qK)*cos(phiS - phiK) + cos(qS)*cos(qK)));
}

Complex polarization(double qS, double phiS, double qK, double phiK){
    // Calculate the rotation of polarization angle :math:`exp(1j*psi)` due to transforming from
    // the plus and cross polarizations in the source frame to the plus and
    // cross polarization in the SSB frame.

    // args:
    //     qS (double): polar angle of the source's sky location
    //     phiS (double): azimuthal angle of the source's sky location
    //     qK (double): polar angle of the Kerr spin vector
    //     phiK (double): azimuthal angle of the Kerr spin vector
    
    double real_part = cos(qS)*sin(qK)*cos(phiK - phiS) - cos(qK)*sin(qS);
    double imag_part = sin(qK)*sin(phiK - phiS);
    if (abs(real_part) + abs(imag_part) == 0.){
		return 0.;
	}

    return pow(Complex(real_part, imag_part), 2)/(pow(real_part, 2) + pow(imag_part, 2));
}

double solar_mass_to_seconds(double mass){
  return mass*GM_const/pow(c_const, 3);
}
double seconds_to_solar_mass(double seconds){
  return seconds/solar_mass_to_seconds(1.);
}
double solar_mass_to_meters(double mass){
  return mass*GM_const/pow(c_const, 2);
}
double solar_mass_to_parsecs(double mass){
  return solar_mass_to_meters(mass)/pc_const;
}
double parsecs_to_solar_mass(double pc){
  return pc/solar_mass_to_parsecs(1.);
}
double seconds_to_years(double seconds){
  return seconds/yr_const;
}
double years_to_seconds(double years){
  return years*yr_const;
}

double scale_strain_amplitude(double mass1, double distance){
  return mass1/parsecs_to_solar_mass(distance*pow(10., 9));
}

double scale_fourier_amplitude(double mass1, double mass2, double distance){
  return solar_mass_to_seconds(mass2)*scale_strain_amplitude(mass1, distance);
}