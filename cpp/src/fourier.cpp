#include "fourier.hpp"

WaveformFourierHarmonicGenerator::WaveformFourierHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts): _Alm(Alm), _mode_selector(Alm, hOpts), _opts(wOpts) {}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformContainer &h, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, HarmonicOptions hOpts, int num_threads, double df, int fsamples){
	HarmonicModeContainer modes = _mode_selector.selectModes(inspiral, theta, hOpts);
	computeWaveformFourierHarmonics(h, modes.lmodes.data(), modes.mmodes.data(), modes.plusY.data(), modes.crossY.data(), modes.lmodes.size(), inspiral, traj, theta, phi, num_threads, df, fsamples);
}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double df, int fsamples){
  double plusY[modeNum];
  double crossY[modeNum];
  for(int i = 0; i < modeNum; i++){
    double sYlm = spin_weighted_spherical_harmonic(-2, l[i], m[i], theta);
    double sYlmMinus = pow(-1, l[i] + m[i])*spin_weighted_spherical_harmonic(2, l[i], m[i], theta);
    plusY[i] = (sYlm + sYlmMinus);
    crossY[i] = (sYlm - sYlmMinus);
  }
  computeWaveformFourierHarmonics(h, l, m, plusY, crossY, modeNum, inspiral, traj, theta, phi, num_threads, df, fsamples);
}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double df, int fsamples){
    double mphi_mod_2pi[modeNum];
    HarmonicSpline2D* Alms[modeNum];
    double twopi = 2.*M_PI;

    // first compute mode-dependent but not time-step dependent information and store
    for(int i = 0; i < modeNum; i++){
      mphi_mod_2pi[i] = fmod(m[i]*phi, twopi);
      Alms[i] = _Alm.getPointer(l[i], m[i]);
    }
    // inspiral will be just over twice as long as the number
    // of frequency samples since we only sample positive
    // frequencies and leave out the zero frequency
    int hmax = h.getSize();
    double dt = inspiral.getTimeSpacing();
	double Tobs = dt*(hmax + 1);
	if(df <= 0.){
      df = 1/Tobs;
    }
    if(fsamples <= 0){
      fsamples = hmax/2;
    }
	
	// std::cout << "Frequency spacing = " << df << "\n";
	// std::cout << "Frequency samples = " << fsamples << "\n";
	// std::cout << "Max frequency = " << df*(fsamples) << "\n";
	// std::cout << "Time spacing = " << dt << "\n";

    double a = inspiral.getSpin();
    double chi = chi_of_spin(a);
    double omega_min = inspiral.getInitialFrequency();
    double omega_max = inspiral.getFinalFrequency();
	double oisco = inspiral.getISCOFrequency();

    double alpha_i = alpha_of_a_omega(a, omega_min);
    double phase_i = traj.phase(chi, alpha_i);
	double time_i = traj.time(chi, alpha_i);
    double massratio = inspiral.getMassRatio();

    // for some reason creating arrays of this size was causing issues with OpenMP
    // minimal internet research suggests that OpenMP allocates arrays on the stack
    // and large arrays can lead to crashes. std::vector seems to work better here
    Vector hplusReal(modeNum*fsamples);
    Vector hplusImag(modeNum*fsamples);
    Vector hcrossReal(modeNum*fsamples);
    Vector hcrossImag(modeNum*fsamples);

    #pragma omp parallel num_threads(num_threads)
    {
      int i, j;
      double amp, modePhase, Phi, cPhi, sPhi, omega, alpha, dtdo, phase, t;
      // first we calculate all of the mode data
      #pragma omp for collapse(2)
      for(j = 0; j < modeNum; j++){
        for(i = 0; i < fsamples; i++){
          omega = twopi/m[j]*(i + 1)*df;
          if(omega_min <= omega && omega <= omega_max){
            alpha = alpha_of_a_omega(a, omega, oisco);
            phase = (traj.phase(chi, alpha) - phase_i)/massratio;
			t = (traj.time(chi, alpha) - time_i)/massratio;
			dtdo = abs(traj.time_of_a_alpha_omega_derivative(a, alpha))/massratio;
            amp = Alms[j]->amplitude(chi, alpha)*sqrt(twopi/m[j]*dtdo);
            modePhase = Alms[j]->phase(chi, alpha);
            Phi = modePhase - fmod(m[j]*phase, twopi) + mphi_mod_2pi[j];
            Phi += fmod(m[j]*omega*t, twopi) - 0.25*M_PI;
            cPhi = std::cos(Phi);
            sPhi = std::sin(Phi);
            hplusReal[j + i*modeNum] = 0.5*amp*plusY[j]*cPhi;
            hplusImag[j + i*modeNum] = 0.5*amp*plusY[j]*sPhi;
            hcrossReal[j + i*modeNum] = -0.5*amp*crossY[j]*sPhi;
            hcrossImag[j + i*modeNum] = 0.5*amp*crossY[j]*cPhi;
          }
        }
      }

      // we then sum over the data and have it stored in the waveform container
      #pragma omp for schedule(static)
      for(int i = 0; i < fsamples; i++){
        for(int j = 0; j < modeNum; j++){
          h.addTimeStep(i, hplusReal[j + i*modeNum], hcrossReal[j + i*modeNum]);
          h.addTimeStep(hmax - 1 - i, hplusImag[j + i*modeNum], hcrossImag[j + i*modeNum]);
        }
      }
    }


}

HarmonicSelector& WaveformFourierHarmonicGenerator::getModeSelector(){
	return _mode_selector;
}

HarmonicModeContainer WaveformFourierHarmonicGenerator::selectModes(InspiralContainer &inspiral, double theta){
	return _mode_selector.selectModes(inspiral, theta);
}

HarmonicModeContainer WaveformFourierHarmonicGenerator::selectModes(InspiralContainer &inspiral, double theta, HarmonicOptions opts){
	return _mode_selector.selectModes(inspiral, theta, opts);
}

WaveformHarmonicOptions WaveformFourierHarmonicGenerator::getWaveformHarmonicOptions(){
	return _opts;
}

HarmonicOptions WaveformFourierHarmonicGenerator::getHarmonicOptions(){
	return _mode_selector.getHarmonicOptions();
}

// Waveform Generator

WaveformFourierGenerator::WaveformFourierGenerator(TrajectorySpline2D &traj, HarmonicAmplitudes &harm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts): WaveformFourierHarmonicGenerator(harm, hOpts, wOpts), _inspiralGen(traj) {}

double WaveformFourierGenerator::convertTime(double t, double M){
  /*
  Converts time in seconds to time in M_\odot
  */
	return t/solar_mass_to_seconds(M);
}

double WaveformFourierGenerator::convertFrequency(double f, double M){
  /*
  Converts frequency in Hz to frequency in 1/M_\odot
  */
	return f*solar_mass_to_seconds(M);
}

int WaveformFourierGenerator::computeFrequencyStepNumber(double dt, double T){
	return years_to_seconds(T)/dt/2;
}

void WaveformFourierGenerator::computeFourierWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_fourier_amplitude(mu, M, dist);

	// omp_set_num_threads(16);
	StopWatch watch;
	// watch.start();
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, wOpts.num_threads);
	// watch.stop();
	// watch.print();
	// watch.reset();
	computeWaveformFourierHarmonics(h, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, hOpts, wOpts.num_threads);
	
	double rescaleRe, rescaleIm;
	rescaleRe = std::real(wOpts.rescale);
	rescaleIm = std::imag(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes
	// watch.start();
	int imax = h.getSize();
  	int imaxf = imax/2;
	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		double hplus, hcross;
		#pragma omp for
		for(i = 0; i < imaxf; i++){
      		// real components of the Fourier harmonics
			hplus = h.getPlus(i);
			hcross = h.getCross(i);
			h.setTimeStep(i, rescaleRe*hplus + rescaleIm*hcross, rescaleRe*hcross - rescaleIm*hplus);
			// imaginary components of the Fourier harmonics
			hplus = h.getPlus(imax - 1 - i);
			hcross = h.getCross(imax - 1 - i);
			h.setTimeStep(imax - 1 - i, rescaleRe*hplus + rescaleIm*hcross, rescaleRe*hcross - rescaleIm*hplus);
		}
	}
	// watch.stop();
	// watch.print();
	// watch.reset();
}

void WaveformFourierGenerator::computeFourierWaveform(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_fourier_amplitude(mu, M, dist);

	// omp_set_num_threads(16);
	StopWatch watch;
	// watch.start();
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, wOpts.num_threads);
	// watch.stop();
	// watch.print();
	// watch.reset();
	computeWaveformFourierHarmonics(h, l, m, modeNum, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, wOpts.num_threads);
	
	double rescaleRe, rescaleIm;
	rescaleRe = std::real(wOpts.rescale);
	rescaleIm = std::imag(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes
	// watch.start();
	int imax = h.getSize();
  	int imaxf = imax/2;
	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		double hplus, hcross;
		#pragma omp for
		for(i = 0; i < imaxf; i++){
			// real components of the Fourier harmonics
			hplus = h.getPlus(i);
			hcross = h.getCross(i);
			h.setTimeStep(i, rescaleRe*hplus + rescaleIm*hcross, rescaleRe*hcross - rescaleIm*hplus);
			// imaginary components of the Fourier harmonics
			hplus = h.getPlus(imax - 1 - i);
			hcross = h.getCross(imax - 1 - i);
			h.setTimeStep(imax - 1 - i, rescaleRe*hplus + rescaleIm*hcross, rescaleRe*hcross - rescaleIm*hplus);
		}
	}
	// watch.stop();
	// watch.print();
	// watch.reset();
}

void WaveformFourierGenerator::computeFourierWaveformSourceFrame(WaveformContainer &h, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T){
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions opts = getWaveformHarmonicOptions();

	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, opts.num_threads);
	computeWaveformFourierHarmonics(h, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, getHarmonicOptions(), opts.num_threads);

	double amplitude_correction = solar_mass_to_seconds(M);
	#pragma omp parallel num_threads(opts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		#pragma omp for
		for(i = 0; i < h.getSize(); i++){
			h.multiplyTimeStep(i, amplitude_correction, amplitude_correction);
		}
	}
}

void WaveformFourierGenerator::computeFourierWaveformSourceFrame(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double dt, double T){
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions opts = getWaveformHarmonicOptions();

	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, opts.num_threads);
	computeWaveformFourierHarmonics(h, l, m, modeNum, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, opts.num_threads);

	double amplitude_correction = solar_mass_to_seconds(M);
	#pragma omp parallel num_threads(opts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		#pragma omp for
		for(i = 0; i < h.getSize(); i++){
			h.multiplyTimeStep(i, amplitude_correction, amplitude_correction);
		}
	}
}

HarmonicModeContainer WaveformFourierGenerator::selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T){
	return selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, dt, T, getHarmonicOptions());
}

HarmonicModeContainer WaveformFourierGenerator::selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double dt, double T, HarmonicOptions opts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions wOpts = getWaveformHarmonicOptions();

	dt = T/(opts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, wOpts.num_threads);
	return WaveformFourierHarmonicGenerator::selectModes(inspiral, theta, opts);
}

double scale_fourier_amplitude(double mass1, double mass2, double distance){
  return solar_mass_to_seconds(mass2)*scale_strain_amplitude(mass1, distance);
}
