#include "fourier.hpp"

WaveformFourierHarmonicGenerator::WaveformFourierHarmonicGenerator(HarmonicAmplitudes &Alm, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts): _Alm(Alm), _mode_selector(Alm, hOpts), _opts(wOpts) {}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformContainer &h, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, HarmonicOptions hOpts, int num_threads, double freq[], int fsamples){
	HarmonicModeContainer modes = _mode_selector.selectModes(inspiral, theta, hOpts);
	computeWaveformFourierHarmonics(h, modes.lmodes.data(), modes.mmodes.data(), modes.plusY.data(), modes.crossY.data(), modes.lmodes.size(), inspiral, traj, theta, phi, num_threads, freq, fsamples);
}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples, int include_negative_m){
  double plusY[modeNum];
  double crossY[modeNum];
  double sYlm, sYlmMinus;
  int mm;
  for(int i = 0; i < modeNum; i++){
	mm = abs(m[i]);
	if(include_negative_m){
		sYlm = spin_weighted_spherical_harmonic(-2, l[i], mm, theta);
    	sYlmMinus = pow(-1, l[i] + mm)*spin_weighted_spherical_harmonic(2, l[i], mm, theta);
    	plusY[i] = (sYlm + sYlmMinus);
    	crossY[i] = (sYlm - sYlmMinus);
	}else if(m[i] > 0){
		sYlm = spin_weighted_spherical_harmonic(-2, l[i], mm, theta);
    	plusY[i] = sYlm;
    	crossY[i] = sYlm;
	}else if(m[i] < 0){
		sYlmMinus = pow(-1, l[i] + mm)*spin_weighted_spherical_harmonic(2, l[i], mm, theta);
    	plusY[i] = sYlmMinus;
    	crossY[i] = -sYlmMinus;
	}else{
		plusY[i] = 0.;
    	crossY[i] = 0.;
	}
  }
  computeWaveformFourierHarmonics(h, l, m, plusY, crossY, modeNum, inspiral, traj, theta, phi, num_threads, freq, fsamples);
}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples, int include_negative_m){
  double plusY[modeNum];
  double crossY[modeNum];
  double sYlm, sYlmMinus;
  int mm;
  for(int i = 0; i < modeNum; i++){
	mm = abs(m[i]);
	if(include_negative_m){
		sYlm = spin_weighted_spherical_harmonic(-2, l[i], mm, theta);
    	sYlmMinus = pow(-1, l[i] + mm)*spin_weighted_spherical_harmonic(2, l[i], mm, theta);
    	plusY[i] = (sYlm + sYlmMinus);
    	crossY[i] = (sYlm - sYlmMinus);
	}else if(m[i] > 0){
		sYlm = spin_weighted_spherical_harmonic(-2, l[i], mm, theta);
    	plusY[i] = sYlm;
    	crossY[i] = sYlm;
	}else if(m[i] < 0){
		sYlmMinus = pow(-1, l[i] + mm)*spin_weighted_spherical_harmonic(2, l[i], mm, theta);
    	plusY[i] = sYlmMinus;
    	crossY[i] = -sYlmMinus;
	}else{
		plusY[i] = 0.;
    	crossY[i] = 0.;
	}
  }
  computeWaveformFourierHarmonics(h, l, m, plusY, crossY, modeNum, inspiral, traj, theta, phi, num_threads, freq, fsamples);
}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples){
    double mphi_mod_2pi[modeNum];
    HarmonicSpline2D* Alms[modeNum];
    double twopi = 2.*M_PI;
	int maxM = 1;
	int minM = 15;

    // first compute mode-dependent but not time-step dependent information and store
    for(int i = 0; i < modeNum; i++){
	  int mm = abs(m[i]);
      mphi_mod_2pi[i] = fmod(mm*phi, twopi);
      Alms[i] = _Alm.getPointer(l[i], mm);
	  if(mm > maxM){maxM = mm;}
	  if(mm < minM){minM = mm;}
    }
    // inspiral will be just over twice as long as the number
    // of frequency samples since we only sample positive
    // frequencies and leave out the zero frequency
    int hmax = h.getSize();
    if(fsamples <= 0){
      fsamples = hmax/2;
    }
	if(fsamples > hmax/2){
		std::cout << "(FOURIER) Error: More frequency samples than waveform samples \n";
	}

    double a = inspiral.getSpin();
    double chi = chi_of_spin(a);
    double omega_min = inspiral.getInitialFrequency();
    double omega_max = inspiral.getFinalFrequency();
	double oisco = inspiral.getISCOFrequency();

	// std::cout << "omega_min = " << omega_min << "\n";
	// std::cout << "omega_max = " << omega_max << "\n";
	// std::cout << "oisco = " << oisco << "\n";
	// std::cout << "alpha_f = " << inspiral.getAlpha(inspiral.getSize() - 1) << "\n";

    double alpha_i = alpha_of_a_omega(a, omega_min);
    double phase_i = traj.phase(chi, alpha_i);
	double time_i = traj.time(chi, alpha_i);
    double massratio = inspiral.getMassRatio();

    // for some reason creating arrays of this size was causing issues with OpenMP
    // minimal internet research suggests that OpenMP allocates arrays on the stack
    // and large arrays can lead to crashes. std::vector seems to work better here
    Vector hplusReal(modeNum*fsamples, 0.);
    Vector hplusImag(modeNum*fsamples, 0.);
    Vector hcrossReal(modeNum*fsamples, 0.);
    Vector hcrossImag(modeNum*fsamples, 0.);


	// attempt to limit some of the for-loops
	// basically assess which frequency samples will be non-zero before trying to
	// calculate them in the series of for-loops below
	double fmin = minM*omega_min/twopi;
	int freq_min_iter = 0;
	while(freq[freq_min_iter] < fmin && freq_min_iter < fsamples - 1){
		freq_min_iter++;
	}

	double fmax = maxM*omega_max/twopi;
	int freq_max_iter = fsamples - 1;
	while(freq[freq_max_iter] > fmax && freq_max_iter > freq_min_iter + 1){
		freq_max_iter--;
	}

	int freq_iter_samples = freq_max_iter - freq_min_iter + 1;

    #pragma omp parallel num_threads(num_threads)
    {
      int i, j, k;
      double amp, modePhase, Phi, cPhi, sPhi, omega, alpha, dtdo, deltaPhase;
      // first we calculate all of the mode data
      #pragma omp for collapse(2)
      for(k = 0; k < freq_iter_samples; k++){
		for(j = 0; j < modeNum; j++){
		  int mm = abs(m[j]);
		  i = freq_min_iter + k;
		  omega = twopi*freq[i]/mm; // from m\omega = 2\pi f
		  if(omega >= omega_min && omega <= omega_max){
			alpha = alpha_of_a_omega(a, omega, oisco);
			deltaPhase = (traj.phase(chi, alpha) - omega*traj.time(chi, alpha) - phase_i + omega*time_i)/massratio;
			dtdo = abs(traj.time_of_a_alpha_omega_derivative(a, alpha))/massratio;
			modePhase = Alms[j]->phase(chi, alpha);
			amp = Alms[j]->amplitude(chi, alpha)*sqrt(twopi/mm*dtdo);

			Phi = modePhase - fmod(mm*deltaPhase, twopi) + mphi_mod_2pi[j] - 0.25*M_PI;
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

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonics(WaveformHarmonicsContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples){
    double mphi_mod_2pi[modeNum];
    HarmonicSpline2D* Alms[modeNum];
    double twopi = 2.*M_PI;
	int maxM = 1;
	int minM = 15;

    // first compute mode-dependent but not time-step dependent information and store
    for(int i = 0; i < modeNum; i++){
	  int mm = abs(m[i]);
      mphi_mod_2pi[i] = fmod(mm*phi, twopi);
      Alms[i] = _Alm.getPointer(l[i], mm);
	  if(mm > maxM){maxM = mm;}
	  if(mm < minM){minM = mm;}
    }
    // inspiral will be just over twice as long as the number
    // of frequency samples since we only sample positive
    // frequencies and leave out the zero frequency
    int hmax = h.getTimeSize();
    if(fsamples <= 0){
      fsamples = hmax/2;
    }
	if(fsamples > hmax/2){
		std::cout << "(FOURIER) Error: More frequency samples than waveform samples \n";
	}

    double a = inspiral.getSpin();
    double chi = chi_of_spin(a);
    double omega_min = inspiral.getInitialFrequency();
    double omega_max = inspiral.getFinalFrequency();
	double oisco = inspiral.getISCOFrequency();

	// std::cout << "omega_min = " << omega_min << "\n";
	// std::cout << "omega_max = " << omega_max << "\n";
	// std::cout << "oisco = " << oisco << "\n";
	// std::cout << "alpha_f = " << inspiral.getAlpha(inspiral.getSize() - 1) << "\n";

    double alpha_i = alpha_of_a_omega(a, omega_min);
    double phase_i = traj.phase(chi, alpha_i);
	double time_i = traj.time(chi, alpha_i);
    double massratio = inspiral.getMassRatio();

    // for some reason creating arrays of this size was causing issues with OpenMP
    // minimal internet research suggests that OpenMP allocates arrays on the stack
    // and large arrays can lead to crashes. std::vector seems to work better here
    Vector hplusReal(modeNum*fsamples, 0.);
    Vector hplusImag(modeNum*fsamples, 0.);
    Vector hcrossReal(modeNum*fsamples, 0.);
    Vector hcrossImag(modeNum*fsamples, 0.);


	// attempt to limit some of the for-loops
	// basically assess which frequency samples will be non-zero before trying to
	// calculate them in the series of for-loops below
	double fmin = minM*omega_min/twopi;
	int freq_min_iter = 0;
	while(freq[freq_min_iter] < fmin && freq_min_iter < fsamples - 1){
		freq_min_iter++;
	}

	double fmax = maxM*omega_max/twopi;
	int freq_max_iter = fsamples - 1;
	while(freq[freq_max_iter] > fmax && freq_max_iter > freq_min_iter + 1){
		freq_max_iter--;
	}

	int freq_iter_samples = freq_max_iter - freq_min_iter + 1;

    #pragma omp parallel num_threads(num_threads)
    {
      int i, j, k;
      double amp, modePhase, Phi, cPhi, sPhi, omega, alpha, dtdo, deltaPhase;
      // first we calculate all of the mode data
      #pragma omp for collapse(2)
      for(k = 0; k < freq_iter_samples; k++){
		for(j = 0; j < modeNum; j++){
		  int mm = abs(m[j]);
		  i = freq_min_iter + k;
		  omega = twopi*freq[i]/mm; // from m\omega = 2\pi f
		  if(omega >= omega_min && omega <= omega_max){
			alpha = alpha_of_a_omega(a, omega, oisco);
			deltaPhase = (traj.phase(chi, alpha) - omega*traj.time(chi, alpha) - phase_i + omega*time_i)/massratio;
			dtdo = abs(traj.time_of_a_alpha_omega_derivative(a, alpha))/massratio;
			modePhase = Alms[j]->phase(chi, alpha);
			amp = Alms[j]->amplitude(chi, alpha)*sqrt(twopi/mm*dtdo);

			Phi = modePhase - fmod(mm*deltaPhase, twopi) + mphi_mod_2pi[j] - 0.25*M_PI;
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
          h.setTimeStep(j, i, hplusReal[j + i*modeNum], hcrossReal[j + i*modeNum]);
          h.setTimeStep(j, hmax - 1 - i, hplusImag[j + i*modeNum], hcrossImag[j + i*modeNum]);
        }
      }
    }
}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonicsPhaseAmplitude(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples){
  double plusY[modeNum];
  double crossY[modeNum];
  double sYlm, sYlmMinus;
  int mm;
  for(int i = 0; i < modeNum; i++){
	mm = abs(m[i]);
	sYlm = spin_weighted_spherical_harmonic(-2, l[i], mm, theta);
	sYlmMinus = pow(-1, mm)*spin_weighted_spherical_harmonic(2, l[i], mm, theta);
	plusY[i] = sYlm;
	crossY[i] = sYlmMinus;
  }
  computeWaveformFourierHarmonicsPhaseAmplitude(h, l, m, plusY, crossY, modeNum, inspiral, traj, theta, phi, num_threads, freq, fsamples);
}

void WaveformFourierHarmonicGenerator::computeWaveformFourierHarmonicsPhaseAmplitude(WaveformHarmonicsContainer &h, int l[], int m[], double plusY[], double crossY[], int modeNum, InspiralContainer &inspiral, TrajectorySpline2D &traj, double theta, double phi, int num_threads, double freq[], int fsamples){
    // double mphi_mod_2pi[modeNum];
    HarmonicSpline2D* Alms[modeNum];
    double twopi = 2.*M_PI;
	int maxM = 1;
	int minM = 15;

    // first compute mode-dependent but not time-step dependent information and store
    for(int i = 0; i < modeNum; i++){
	  int mm = abs(m[i]);
    //   mphi_mod_2pi[i] = fmod(mm*phi, twopi);
      Alms[i] = _Alm.getPointer(l[i], mm);
	  if(mm > maxM){maxM = mm;}
	  if(mm < minM){minM = mm;}
    }
    // inspiral will be just over twice as long as the number
    // of frequency samples since we only sample positive
    // frequencies and leave out the zero frequency
    int hmax = h.getTimeSize();
    if(fsamples <= 0){
      fsamples = hmax;
    }
	if(fsamples > hmax){
		std::cout << "(FOURIER) Error: More frequency samples than waveform samples \n";
	}

    double a = inspiral.getSpin();
    double chi = chi_of_spin(a);
    double omega_min = inspiral.getInitialFrequency();
    double omega_max = inspiral.getFinalFrequency();
	double oisco = inspiral.getISCOFrequency();

	// std::cout << "omega_min = " << omega_min << "\n";
	// std::cout << "omega_max = " << omega_max << "\n";
	// std::cout << "oisco = " << oisco << "\n";
	// std::cout << "alpha_f = " << inspiral.getAlpha(inspiral.getSize() - 1) << "\n";

    double alpha_i = alpha_of_a_omega(a, omega_min);
    double phase_i = traj.phase(chi, alpha_i);
	double time_i = traj.time(chi, alpha_i);
    double massratio = inspiral.getMassRatio();

    // for some reason creating arrays of this size was causing issues with OpenMP
    // minimal internet research suggests that OpenMP allocates arrays on the stack
    // and large arrays can lead to crashes. std::vector seems to work better here
    // Vector hplusAmp(modeNum*fsamples, 0.);
    // Vector hplusPhase(modeNum*fsamples, 0.);
    // Vector hminusAmp(modeNum*fsamples, 0.);
    // Vector hminusPhase(modeNum*fsamples, 0.);


	// attempt to limit some of the for-loops
	// basically assess which frequency samples will be non-zero before trying to
	// calculate them in the series of for-loops below
	double fmin = minM*omega_min/twopi;
	int freq_min_iter = 0;
	while(freq[freq_min_iter] < fmin && freq_min_iter < fsamples - 1){
		freq_min_iter++;
	}

	double fmax = maxM*omega_max/twopi;
	int freq_max_iter = fsamples - 1;
	while(freq[freq_max_iter] > fmax && freq_max_iter > freq_min_iter + 1){
		freq_max_iter--;
	}

	int freq_iter_samples = freq_max_iter - freq_min_iter + 1;

    #pragma omp parallel num_threads(num_threads)
    {
      int i, j, k;
      double amp, modePhase, Phi, cPhi, sPhi, omega, alpha, dtdo, deltaPhase;
      // first we calculate all of the mode data
      #pragma omp for collapse(2)
      for(k = 0; k < freq_iter_samples; k++){
		for(j = 0; j < modeNum; j++){
		  int mm = abs(m[j]);
		  i = freq_min_iter + k;
		  omega = twopi*freq[i]/mm; // from m\omega = 2\pi f
		  if(omega >= omega_min && omega <= omega_max){
			alpha = alpha_of_a_omega(a, omega, oisco);
			deltaPhase = (traj.phase(chi, alpha) - omega*traj.time(chi, alpha) - phase_i + omega*time_i)/massratio;
			dtdo = abs(traj.time_of_a_alpha_omega_derivative(a, alpha))/massratio;
			modePhase = Alms[j]->phase(chi, alpha);
			amp = Alms[j]->amplitude(chi, alpha)*sqrt(twopi/mm*dtdo);
			Phi = modePhase - mm*(deltaPhase - phi) - 0.25*M_PI;
			// Phi = modePhase - fmod(mm*deltaPhase, twopi) + mphi_mod_2pi[j] - 0.25*M_PI;

			h.setTimeStep(2*j, i, amp*plusY[j], Phi); // \tilde{h}_{lm}(f)
			h.setTimeStep(2*j + 1, i, amp*crossY[j],  - Phi - l[j]*M_PI); // \tilde{h}_{l-m}(-f)
			// hplusAmp[j + i*modeNum] = 0.5*amp*plusY[j];
			// hplusPhase[j + i*modeNum] = Phi;
			// hminusAmp[j + i*modeNum] = 0.5*amp*crossY[j];
			// hminusPhase[j + i*modeNum] = - Phi - l[j]*M_PI;
		  }
        }
      }

      // we then sum over the data and have it stored in the waveform container
    //   #pragma omp for schedule(static)
    //   for(int i = 0; i < fsamples; i++){
    //     for(int j = 0; j < modeNum; j++){
    //       h.setTimeStep(j, i, hplusAmp[j + i*modeNum], hplusPhase[j + i*modeNum]);
    //       h.setTimeStep(j, hmax - 1 - i, hminusAmp[j + i*modeNum], hminusPhase[j + i*modeNum]);
    //     }
    //   }
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
	return (years_to_seconds(T)/dt + 1)/2; // we do not add one because we leave off the zero mode
}

int WaveformFourierGenerator::computeTimeStepNumber(double dt, double T){
	return years_to_seconds(T)/dt + 1;
}

int WaveformFourierGenerator::computeTimeStepNumber(double M, double mu, double a, double r0, double dt, double T){
	dt = convertTime(dt, M);
	T = convertTime(years_to_seconds(T), M);
	return _inspiralGen.computeTimeStepNumber(a, mu/M, r0, dt, T);
}

void WaveformFourierGenerator::computeFourierWaveform(WaveformContainer &h, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	
	int imax = h.getSize();
  	int imaxf = imax/2;
	double df;
	Vector freq(imaxf, 0.);

	for(int i = 0; i < imaxf; i ++){
		df = frequencies[i];
		freq[i] = convertFrequency(df, M);
	}
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_fourier_amplitude(mu, M, dist);

	// omp_set_num_threads(16);
	// StopWatch watch;
	// watch.start();
	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(hOpts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	// watch.stop();
	// watch.print();
	// watch.reset();
	computeWaveformFourierHarmonics(h, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, hOpts, wOpts.num_threads, &freq[0], imaxf);
	
	double rescaleRe, rescaleIm;
	rescaleRe = std::real(wOpts.rescale);
	rescaleIm = std::imag(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes
	// watch.start();
	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		double hplusRe, hcrossRe, hplusIm, hcrossIm;
		#pragma omp for
		for(i = 0; i < imaxf; i++){
      		// real components of the Fourier harmonics
			hplusRe = h.getPlus(i);
			hcrossRe = h.getCross(i);
			hplusIm = h.getPlus(imax - 1 - i);
			hcrossIm = h.getCross(imax - 1 - i);
			h.setTimeStep(i, rescaleRe*hplusRe + rescaleIm*hcrossRe, rescaleRe*hcrossRe - rescaleIm*hplusRe);
			h.setTimeStep(imax - 1 - i, rescaleRe*hplusIm + rescaleIm*hcrossIm, rescaleRe*hcrossIm - rescaleIm*hplusIm);
			// h.setTimeStep(i, rescaleRe*hplusRe - rescaleIm*hplusIm, rescaleRe*hcrossRe - rescaleIm*hcrossIm);
			// h.setTimeStep(imax - 1 - i, rescaleRe*hplusIm + rescaleIm*hplusRe, rescaleRe*hcrossIm + rescaleIm*hcrossRe);
		}
	}
	// watch.stop();
	// watch.print();
	// watch.reset();
}

void WaveformFourierGenerator::computeFourierWaveform(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);

	int imax = h.getSize();
  	int imaxf = imax/2;
	double fi;
	Vector freq(imaxf, 0.);

	for(int i = 0; i < imaxf; i ++){
		fi = frequencies[i];
		freq[i] = convertFrequency(fi, M);
	}
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_fourier_amplitude(mu, M, dist);

	// omp_set_num_threads(16);
	// StopWatch watch;
	// watch.start();
	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(hOpts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	// watch.stop();
	// watch.print();
	// watch.reset();
	computeWaveformFourierHarmonics(h, l, m, modeNum, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, wOpts.num_threads, &freq[0], imaxf, wOpts.include_negative_m);
	
	double rescaleRe, rescaleIm;
	rescaleRe = std::real(wOpts.rescale);
	rescaleIm = std::imag(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes
	// watch.start();

	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		// total_td = omp_get_num_threads();
		int i;
		double hplusRe, hcrossRe, hplusIm, hcrossIm;
		#pragma omp for
		for(i = 0; i < imaxf; i++){
      		// real components of the Fourier harmonics
			hplusRe = h.getPlus(i);
			hcrossRe = h.getCross(i);
			hplusIm = h.getPlus(imax - 1 - i);
			hcrossIm = h.getCross(imax - 1 - i);
			h.setTimeStep(i, rescaleRe*hplusRe + rescaleIm*hcrossRe, rescaleRe*hcrossRe - rescaleIm*hplusRe);
			h.setTimeStep(imax - 1 - i, rescaleRe*hplusIm + rescaleIm*hcrossIm, rescaleRe*hcrossIm - rescaleIm*hplusIm);
			// h.setTimeStep(i, rescaleRe*hplusRe - rescaleIm*hplusIm, rescaleRe*hcrossRe - rescaleIm*hcrossIm);
			// h.setTimeStep(imax - 1 - i, rescaleRe*hplusIm + rescaleIm*hplusRe, rescaleRe*hcrossIm + rescaleIm*hcrossRe);
		}
	}
	// watch.stop();
	// watch.print();
	// watch.reset();
}

void WaveformFourierGenerator::computeFourierWaveform(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);

	int imax = h.getTimeSize();
  	int imaxf = imax/2;
	double fi;
	Vector freq(imaxf, 0.);

	for(int i = 0; i < imaxf; i ++){
		fi = frequencies[i];
		freq[i] = convertFrequency(fi, M);
	}
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_fourier_amplitude(mu, M, dist);

	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(hOpts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	computeWaveformFourierHarmonics(h, l, m, modeNum, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, wOpts.num_threads, &freq[0], imaxf, wOpts.include_negative_m);
	
	double rescaleRe, rescaleIm;
	rescaleRe = std::real(wOpts.rescale);
	rescaleIm = std::imag(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes

	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		int i;
		double hplusRe, hcrossRe, hplusIm, hcrossIm;
		#pragma omp for
		for(i = 0; i < imaxf; i++){
			for(int j = 0; j < modeNum; j++){
				hplusRe = h.getPlus(j, i);
				hcrossRe = h.getCross(j, i);
				hplusIm = h.getPlus(j, imax - 1 - i);
				hcrossIm = h.getCross(j, imax - 1 - i);
				h.setTimeStep(j, i, rescaleRe*hplusRe + rescaleIm*hcrossRe, rescaleRe*hcrossRe - rescaleIm*hplusRe);
				h.setTimeStep(j, imax - 1 - i, rescaleRe*hplusIm + rescaleIm*hcrossIm, rescaleRe*hcrossIm - rescaleIm*hplusIm);
			}
		}
	}
}

void WaveformFourierGenerator::computeFourierWaveformPhaseAmplitude(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double dist, double qS, double phiS, double qK, double phiK, double Phi_phi0, double frequencies[], double T, HarmonicOptions hOpts, WaveformHarmonicOptions wOpts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);

	int imax = h.getTimeSize();
	int jmax = h.getModeSize();
	// if(jmax != 2*modeNum){
	// 	std::cout << "ERROR";
	// }
  	int imaxf = imax;
	double fi;
	Vector freq(imaxf, 0.);

	for(int i = 0; i < imaxf; i ++){
		fi = frequencies[i];
		freq[i] = convertFrequency(fi, M);
	}
	T = convertTime(years_to_seconds(T), M);
	wOpts.rescale = polarization(qS, phiS, qK, phiK);
	wOpts.rescale *= scale_fourier_amplitude(mu, M, dist);

	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(hOpts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	computeWaveformFourierHarmonicsPhaseAmplitude(h, l, m, modeNum, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, wOpts.num_threads, &freq[0], imaxf);
	
	double rescaleAmp, rescalePhase;
	rescaleAmp = std::abs(wOpts.rescale);
	rescalePhase = std::arg(wOpts.rescale);
	// if the rescaling factor is purely real, then just rescale both polarizations by the same amplitude
	// else, then we get a mixing of the plus and cross polarizations that gives us new polarization amplitudes

	#pragma omp parallel num_threads(wOpts.num_threads)
	{
		int i, j;
		double amp, phase;
		#pragma omp for
		for(i = 0; i < imaxf; i++){
			for(j = 0; j < jmax; j++){
				amp = h.getPlus(j, i);
				if(fabs(amp) > 0.){
					phase = h.getCross(j, i);
					h.setTimeStep(j, i, rescaleAmp*amp, phase + rescalePhase);
				}
			}
		}
	}
}

void WaveformFourierGenerator::computeFourierWaveformSourceFrame(WaveformContainer &h, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double frequencies[], double T){
	int imax = h.getSize();
  	int imaxf = imax/2;
	double df;
	Vector freq(imaxf, 0.);

	for(int i = 0; i < imaxf; i ++){
		df = frequencies[i];
		freq[i] = convertFrequency(df, M);
	}
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions opts = getWaveformHarmonicOptions();
	HarmonicOptions hOpts = getHarmonicOptions();

	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(hOpts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	computeWaveformFourierHarmonics(h, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, getHarmonicOptions(), opts.num_threads, &freq[0], imaxf);

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

void WaveformFourierGenerator::computeFourierWaveformSourceFrame(WaveformContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double frequencies[], double T){
	int imax = h.getSize();
  	int imaxf = imax/2;
	double df;
	Vector freq(imaxf, 0.);

	for(int i = 0; i < imaxf; i ++){
		df = frequencies[i];
		freq[i] = convertFrequency(df, M);
	}
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions opts = getWaveformHarmonicOptions();
	HarmonicOptions hOpts = getHarmonicOptions();

	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(hOpts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	computeWaveformFourierHarmonics(h, l, m, modeNum, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, opts.num_threads, &freq[0], imaxf, opts.include_negative_m);

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

void WaveformFourierGenerator::computeFourierWaveformSourceFrame(WaveformHarmonicsContainer &h, int l[], int m[], int modeNum, double M, double mu, double a, double r0, double theta, double phi, double Phi_phi0, double frequencies[], double T){
	int imax = h.getTimeSize();
  	int imaxf = imax/2;
	double df;
	Vector freq(imaxf, 0.);

	for(int i = 0; i < imaxf; i ++){
		df = frequencies[i];
		freq[i] = convertFrequency(df, M);
	}
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions opts = getWaveformHarmonicOptions();
	HarmonicOptions hOpts = getHarmonicOptions();

	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(hOpts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	computeWaveformFourierHarmonics(h, l, m, modeNum, inspiral, _inspiralGen.getTrajectorySpline(), theta, phi - Phi_phi0, opts.num_threads, &freq[0], imaxf, opts.include_negative_m);

	double amplitude_correction = solar_mass_to_seconds(M);
	#pragma omp parallel num_threads(opts.num_threads)
	{
		int i, j;
		#pragma omp for
		for(i = 0; i < imax; i++){
			for(j = 0; j < modeNum; j++){
				h.multiplyTimeStep(j, i, amplitude_correction, amplitude_correction);
			}
		}
	}
}

HarmonicModeContainer WaveformFourierGenerator::selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double T){
	return selectModes(M, mu, a, r0, qS, phiS, qK, phiK, Phi_phi0, T, getHarmonicOptions());
}

HarmonicModeContainer WaveformFourierGenerator::selectModes(double M, double mu, double a, double r0, double qS, double phiS, double qK, double phiK, double Phi_phi0, double T, HarmonicOptions opts){
	double theta, phi;
	sourceAngles(theta, phi, qS, phiS, qK, phiK);
	T = convertTime(years_to_seconds(T), M);
	WaveformHarmonicOptions wOpts = getWaveformHarmonicOptions();

	double Tmerge = _inspiralGen.computeTimeToMerger(a, mu/M, r0);
	T = (T > Tmerge) ? Tmerge : T;
	double dt = T/(opts.max_samples - 1);
	InspiralContainer inspiral = _inspiralGen.computeInspiral(a, mu/M, r0, dt, T, 1);
	return WaveformFourierHarmonicGenerator::selectModes(inspiral, theta, opts);
}

double scale_fourier_amplitude(double mass1, double mass2, double distance){
  return solar_mass_to_seconds(mass2)*scale_strain_amplitude(mass1, distance);
}
