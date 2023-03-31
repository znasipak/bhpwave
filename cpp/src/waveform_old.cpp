#include "waveform_old.hpp"

#define REL_MODE_EPS 5.e-4
#define OMEGA_DIFF_EPS 1.e-5
#define ALPHA_MAX 1.
#define ALPHA_MIN 0.

TrajectorySpline2D _traj_(read_trajectory_data());

// std::vector<HarmonicSpline2D> _modes_ {
//   HarmonicSpline2D(read_harmonic_mode_data(2, 1)), 
//   HarmonicSpline2D(read_harmonic_mode_data(2, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(3, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(3, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(3, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(4, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(4, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(4, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(4, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(5, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(5, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(5, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(5, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(5, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(6, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(6, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(6, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(6, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(6, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(6, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(7, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(7, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(7, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(7, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(7, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(7, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(7, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(8, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(9, 9)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 9)),
//   HarmonicSpline2D(read_harmonic_mode_data(10, 10)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 9)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 10)),
//   HarmonicSpline2D(read_harmonic_mode_data(11, 11)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 9)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 10)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 11)),
//   HarmonicSpline2D(read_harmonic_mode_data(12, 12)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 9)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 10)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 11)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 12)),
//   HarmonicSpline2D(read_harmonic_mode_data(13, 13)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 9)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 10)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 11)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 12)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 13)),
//   HarmonicSpline2D(read_harmonic_mode_data(14, 14)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 1)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 2)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 3)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 4)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 5)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 6)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 7)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 8)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 9)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 10)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 11)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 12)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 13)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 14)),
//   HarmonicSpline2D(read_harmonic_mode_data(15, 15))
// };

std::vector<HarmonicSpline2D> _modes_ {
  HarmonicSpline2D(read_harmonic_mode_data(2, 1)), 
  HarmonicSpline2D(read_harmonic_mode_data(2, 2))
};
HarmonicSpline2D mode_lm(int l, int m){
  int lmIter = (l*l - l - 2)/2 + m - 1;
  return _modes_[lmIter];
}

// void load_waveform_data(){
//   // std::cout << "Loading data \n";
//   for(int l = 3; l <= 15; l++){
//     for(int m = 1; m <= l; m++){
//       _modes_.push_back(HarmonicSpline2D(read_harmonic_mode_data(l, m)));
//     }
//   }
//   // std::cout << "Mode data loaded \n";
// }

// for(int l = 3; l <= 15; l++){
//   for(int m = 1; m <= l; m++){
//     _modes_.push_back(HarmonicSpline2D(read_harmonic_mode_data(l, m)));
//   }
// }

GWStrain::GWStrain(int N): plus(N, 0.), cross(N, 0.) {}
GWStrainFourier::GWStrainFourier(int N): plusR(N, 0.), plusI(N, 0.), crossR(N, 0.), crossI(N, 0.) {}

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

double time_spacing(double mass1, double mass2, double spin, double r0, double duration_yr){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;

  double chi = chi_of_spin(spin);
  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i);
  double t_i = _traj_.time(chi, alpha_i);
  if(duration > -t_i){
    duration = -t_i;
  }
  double t_f = t_i + duration - 1.e-15;
  double alpha_f = _traj_.orbital_alpha(chi, t_f);
  double omega_f = _traj_.orbital_frequency(spin, t_f);
  double f_nyquist = omega_f/solar_mass_to_seconds(mass2)/(2.*M_PI);

  int lmaxmode = 15;
  List maxmodes(lmaxmode);
  harmonic_selection(maxmodes, chi, alpha_i, alpha_f, REL_MODE_EPS);
  int i = 0;
  while(maxmodes[i] > 0 && i < lmaxmode){
    i++;
  }
  int max_m = i;

  return 0.5/f_nyquist/max_m;
}

double frequency_spacing(double mass1, double mass2, double spin, double r0, double &duration_yr){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double period = solar_mass_to_seconds(duration/massratio*mass2);
  duration_yr = seconds_to_years(period);

  return pow(10, 9)/period;
}

double scale_strain_amplitude(double mass1, double distance){
  return mass1/parsecs_to_solar_mass(distance*pow(10., 9));
}

double scale_fourier_amplitude(double mass1, double mass2, double distance){
  return solar_mass_to_seconds(mass2)*scale_strain_amplitude(mass1, distance);
}

////////////////////////////
// Time Domain Trajectory //
////////////////////////////

double time_to_inspiral(double mass1, double mass2, double spin, double r0){
  double chi = chi_of_spin(spin);
  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i, _traj_.orbital_frequency_isco(chi));
  double t_i = -_traj_.time(chi, alpha_i);

  return mass2*solar_mass_to_seconds(mass2)*t_i/mass1;
}

void generate_trajectory_td(double mass1, double mass2, double spin, double r0, double phi0, double duration_yr, double deltaT_sec_temp){
  double deltaT_sec = deltaT_sec_temp;
  if(deltaT_sec == 0.){
    deltaT_sec = 0.9*time_spacing(mass1, mass2, spin, r0, duration_yr);
  }

  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaT = massratio*seconds_to_solar_mass(deltaT_sec)/mass2;

  // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  double chi = chi_of_spin(spin);
  double omega_ISCO = _traj_.orbital_frequency_isco(chi);
  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i, omega_ISCO);
  double t_i = _traj_.time(chi, alpha_i);
  if(duration > -t_i){
    duration = -t_i;
  }

  int time_steps = duration/deltaT + 1;
  Vector rp(time_steps);
  Vector phase(time_steps);
  double phase_i = _traj_.phase(chi, _traj_.orbital_alpha(chi, t_i));
  for(int j = 0; j < time_steps; j++){
    double alpha = _traj_.orbital_alpha(chi, t_i + deltaT*j);
    double omega = omega_of_a_alpha(spin, alpha, omega_ISCO);
    rp[j] = kerr_geo_radius_circ(spin, omega);
    phase[j] = (_traj_.phase(chi, alpha) - phase_i)/massratio + phi0;
  }

  std::cout << "Generating trajectory with parameters:\n m = "<<mass1<<" M_odot\n M = "<<mass2<<" M_odot\n r0 = "<<r0<<" M\n phi0 = "<<phi0<<"\n T = "<<duration_yr<<" yrs\n dt = "<<deltaT_sec<<" sec\n";
  char buff[100];
  if(deltaT_sec_temp == 0.){
    sprintf(buff, "waveform/x_logm%.4f_logM%.4f_a%.4f_r%.4f_ph%.4f_T%.3f.txt", log10(mass1), log10(mass2), spin, r0, phi0, duration_yr);
  }else{
    sprintf(buff, "waveform/x_logm%.4f_logM%.4f_a%.4f_r%.4f_ph%.4f_T%.3f_dt%.1f.txt", log10(mass1), log10(mass2), spin, r0, phi0, duration_yr, deltaT_sec);
  }
  std::string filepath = buff;

  std::cout << "Exporting trajectory to " << filepath << " \n";

  std::ofstream file;
  file.open(filepath);
  file << "t\tr0\tphi\n";
  for(size_t i = 0; i < rp.size(); i++){
    sprintf(buff, "%.8e\t%.8e\t%.8e\n", i*deltaT_sec, rp[i], phase[i]);
    file << buff;
  }

  file.close();

}

void output_flux(double *flux, double *a, double *omega, int sampleN){
  for(int i = 0; i < sampleN; i++){
    // double omega = kerr_geo_time_frequency_circ(a[i], r[i]);
    flux[i] = _traj_.flux_of_a_omega(a[i], omega[i]);
  }
}

void output_trajectory_td(double *rp, double *phase, double mass1, double mass2, double spin, double r0, double phi0, double duration_yr, double deltaT_sec_temp){
  double deltaT_sec = deltaT_sec_temp;
  if(deltaT_sec == 0.){
    deltaT_sec = 0.9*time_spacing(mass1, mass2, spin, r0, duration_yr);
  }

  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaT = massratio*seconds_to_solar_mass(deltaT_sec)/mass2;

  // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  double chi = chi_of_spin(spin);
  double omega_ISCO = _traj_.orbital_frequency_isco(chi);
  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i, omega_ISCO);
  double t_i = _traj_.time(chi, alpha_i);
  if(duration > -t_i){
    duration = -t_i;
  }

  int time_steps = duration/deltaT + 1;
  double phase_i = _traj_.phase(chi, _traj_.orbital_alpha(chi, t_i));
  for(int j = 0; j < time_steps; j++){
    double alpha = _traj_.orbital_alpha(chi, t_i + deltaT*j);
    double omega = omega_of_a_alpha(spin, alpha, omega_ISCO);
    rp[j] = kerr_geo_radius_circ(spin, omega);
    phase[j] = (_traj_.phase(chi, alpha) - phase_i)/massratio + phi0;
  }

}

///////////////////////////
// Time Domain Waveforms //
///////////////////////////

// void reduced_waveform_harmonic_td(FloatVector &hp, FloatVector &hc, const double &chi, const Vector &alpha, const int &L, const int &m, const double &theta, const double &phi, const Vector &phase){
//   double sYlm = spin_weighted_spherical_harmonic(-2, L, m, theta);
//   double sYlmMinus = spin_weighted_spherical_harmonic(2, L, m, theta);
//   // HarmonicMode modeData = read_harmonic_mode_data(L, m);
//   double mphi_mod_2pi = fmod(m*phi, 2.*M_PI);

//   float amp = mode_lm(L, m).amplitude(chi, alpha[0]);
//   float modePhase = mode_lm(L, m).phase(chi, alpha[0]);
//   float Phi = modePhase - fmod(m*phase[0], 2.*M_PI) + mphi_mod_2pi;
//   hp[0] += amp*(sYlm + pow(-1, L + m)*sYlmMinus)*cos(Phi);
//   hc[0] += -amp*(sYlm - pow(-1, L + m)*sYlmMinus)*sin(Phi);

//   for(size_t i = 1; i < hp.size(); i++){
//     amp = mode_lm(L, m).amplitude(chi, alpha[i]);
//     modePhase = mode_lm(L, m).phase(chi, alpha[i]);
//     Phi = modePhase - fmod(m*phase[i], 2.*M_PI) + mphi_mod_2pi;
//     hp[i] += amp*(sYlm + pow(-1, L + m)*sYlmMinus)*cos(Phi);
//     hc[i] += -amp*(sYlm - pow(-1, L + m)*sYlmMinus)*sin(Phi);
//   }
// }

void reduced_waveform_harmonic_td(FloatVector &hp, FloatVector &hc, double &chi, Vector &alpha, int &L, int &m, double &theta, double &phi, Vector &phase, int num_threads){
  double sYlm = spin_weighted_spherical_harmonic(-2, L, m, theta);
  double sYlmMinus = spin_weighted_spherical_harmonic(2, L, m, theta);
  // HarmonicMode modeData = read_harmonic_mode_data(L, m);
  double mphi_mod_2pi = fmod(m*phi, 2.*M_PI);
  HarmonicSpline2D modeLm = mode_lm(L, m);

  double amp = modeLm.amplitude(chi, alpha[0]);
  double modePhase = modeLm.phase(chi, alpha[0]);
  double Phi = modePhase - fmod(m*phase[0], 2.*M_PI) + mphi_mod_2pi;
  Complex hphase = exp(Complex(0., -Phi));
  hp[0] += amp*(sYlm + pow(-1, L + m)*sYlmMinus)*std::real(hphase);
  hc[0] += amp*(sYlm - pow(-1, L + m)*sYlmMinus)*std::imag(hphase);
  int imax = hp.size();

  if(num_threads > 0){
    omp_set_num_threads(num_threads);
  }
  #pragma omp parallel shared(hp, hc, sYlm, sYlmMinus, mphi_mod_2pi, modeLm, chi, alpha, phase) private(amp, modePhase, Phi, hphase)
  {
    #pragma omp for
      for(int i = 1; i < imax; i++){
        amp = modeLm.amplitude(chi, alpha[i]);
        modePhase = modeLm.phase(chi, alpha[i]);
        Phi = modePhase - fmod(m*phase[i], 2.*M_PI) + mphi_mod_2pi;
        hphase = exp(Complex(0., -Phi));
        hp[i] += amp*(sYlm + pow(-1, L + m)*sYlmMinus)*std::real(hphase);
        hc[i] += amp*(sYlm - pow(-1, L + m)*sYlmMinus)*std::imag(hphase);
      }
  }
}

void reduced_waveform_harmonic_td(FloatVector &hp, FloatVector &hc, double &chi, Vector &alpha, List &maxmodes, double &theta, double &phi, Vector &phase, int num_threads){
  int lmaxmode = maxmodes.size();
  int i = 0;
  while(maxmodes[i] > 0 && i < lmaxmode){
    i++;
  }
  int mmax = i;
  int imax = hp.size();

  for(int m = 1; m <= mmax; m++){
    double mphi_mod_2pi = fmod(m*phi, 2.*M_PI);
    int lmin = (m < 2) ? 2 : m;
    std::vector<float> sYlmSym(maxmodes[m-1]+1);
    std::vector<float> sYlmAntiSym(maxmodes[m-1]+1);
    // std::vector<HarmonicSpline2D> modeLm(maxmodes[m-1]+1, mode_lm(2,2));
    // construct angular dependence outside of the parallelized loop because this should be quick
    for(int l = lmin; l <= maxmodes[m-1]; l++){
      float sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
      float sYlmMinus = pow(-1, l + m)*spin_weighted_spherical_harmonic(2, l, m, theta);
      sYlmSym[l] = sYlm + sYlmMinus;
      sYlmAntiSym[l] = sYlm - sYlmMinus;
      // modeLm[l] = mode_lm(l, m);
    }
    if(num_threads > 0){
      omp_set_num_threads(num_threads);
    }
    #pragma omp parallel shared(imax, maxmodes, hp, hc, chi, alpha, phase, mmax)
    {
      #pragma omp for collapse(2)
        // iterate over l-values for fixed m
        for(int l = lmin; l <= maxmodes[m-1]; l++){
          // iterate over time steps
          for(int i = 0; i < imax; i++){
            float amp = mode_lm(l, m).amplitude(chi, alpha[i]);
            float modePhase = mode_lm(l, m).phase(chi, alpha[i]);
            // float amp = 1.;
            // float modePhase = 0.;
            float Phi = modePhase - fmod(m*phase[i], 2.*M_PI) + mphi_mod_2pi; 
            // Complex hphase = exp(Complex(0., -Phi));
            hp[i] += amp*sYlmSym[l]*std::cos(Phi);
            hc[i] += -amp*sYlmAntiSym[l]*std::sin(Phi);
          }
        }
    }
  }
}

GWStrain waveform_harmonic_td(int L, int m, double theta, double phi, double mass1, double mass2, double spin, double r0, double duration_yr, double deltaT_sec, int num_threads){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaT = massratio*seconds_to_solar_mass(deltaT_sec)/mass2;
  double chi = chi_of_spin(spin);

  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i, _traj_.orbital_frequency_isco(chi));
  double t_i = _traj_.time(chi, alpha_i);
  if(duration > -t_i){
    duration = -t_i;
  }

  int time_steps = duration/deltaT + 1;
  Vector alpha(time_steps);
  Vector phase(time_steps);
  double phase_i = _traj_.phase(chi, _traj_.orbital_alpha(chi, t_i));
  if(num_threads > 0){
    omp_set_num_threads(num_threads);
  }
  #pragma omp parallel shared(_traj_, chi, alpha, phase, time_steps, massratio, deltaT, phase_i)
  {
    #pragma omp for
      for(int j = 0; j < time_steps; j++){
        alpha[j] = _traj_.orbital_alpha(chi, t_i + deltaT*j);
        phase[j] = (_traj_.phase(chi, alpha[j]) - phase_i)/massratio;
    }
  }

  GWStrain h(time_steps);

  reduced_waveform_harmonic_td(h.plus, h.cross, chi, alpha, L, m, theta, phi, phase, num_threads);

  return h;
}

double relative_harmonic_power(int l, int m, double chi, double alphaInitial, double deltaAlpha, Vector &alphaDot){
  double power = 0.;
  int steps = alphaDot.size();
  for(int i = 0; i < steps; i++){
    double alpha = alphaInitial + i*deltaAlpha;
    if(alpha > ALPHA_MAX){ alpha = ALPHA_MAX; }
    if(alpha > ALPHA_MIN){ power += abs(pow(mode_lm(l, m).amplitude(chi, alpha), 2)/alphaDot[i]*deltaAlpha); }
  }
  return power;
}

void harmonic_selection(List &maxmodes, double chi, double alphaInitial, double alphaFinal, double epsilon){
  int steps = 50;
  double deltaAlpha = (alphaFinal - alphaInitial)/(steps - 1);
  double spin = spin_of_chi(chi);
  // double oISCO = _traj_.orbital_frequency_isco(chi);
  // double beta = oISCO*(1. - spin*oISCO);
  Vector alphaDot(steps);
  for(int i = 0; i < steps; i++){
    double alpha = alphaInitial + i*deltaAlpha;
    if(alpha > ALPHA_MAX){ alpha = ALPHA_MAX; }
    if(alpha < ALPHA_MIN){ alpha = ALPHA_MIN; }
    alphaDot[i] = abs(_traj_.orbital_alpha_derivative(chi, _traj_.time(chi, alpha)));
  }

  double power22 = relative_harmonic_power(2, 2, chi, alphaInitial, deltaAlpha, alphaDot);

  maxmodes[0] = 2; // initialize m = 1 to include 1 <= l <= lmax = 2
  maxmodes[1] = 2; // initialize m = 2 to include 2 <= l <= lmax = 2
  int lmodemax = maxmodes.size();
  for(int m = 1; m <= 2; m++){
    for(int l = 3; l <= lmodemax; l++){
      double power = relative_harmonic_power(l, m, chi, alphaInitial, deltaAlpha, alphaDot);
      if(power/power22 > epsilon){
        // if there is significant power in the (l, m)-mode extend lmax of the m-mode to l
        maxmodes[m-1] = l; 
      }else{
        l = lmodemax + 1;
      }
    }
  }

  for(int m = 3; m <= lmodemax; m++){
    for(int l = m; l <= lmodemax; l++){
      double power = relative_harmonic_power(l, m, chi, alphaInitial, deltaAlpha, alphaDot);
      if(power/power22 > epsilon){
        maxmodes[m-1] = l;
      }else{
        if(l == m){
          // if there is not enough power in the l = m mode, then we don't expect any more power in higher modes
          // so we just end the mode selection process
          m = lmodemax + 1;
          // Additionally we remove this m-mode from the summation
          maxmodes[m-1] = 0; 
        }
        l = lmodemax + 1;
      }
    }
  }

  // for(int l = 3; l <= 15; l++){
  //   for(int m = l; m > 0; m--){
  //     double power = relative_harmonic_power(l, m, chi, alphaInitial, deltaAlpha, alphaDot);
  //     // std::cout << "Power for ("<<l<<","<<m<<")-mode = "<<power<<" \n";
  //     if(power/power22 > epsilon){
  //       lmodes.push_back(l);
  //       mmodes.push_back(m);
  //     }else{
  //       if(l == m){
  //         l = 15;
  //       }
  //       m = 0;
  //     }
  //   }
  // }
}

int inspiral_time_steps(double mass1, double mass2, double spin, double r0, double duration_yr, double deltaT_sec){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaT = massratio*seconds_to_solar_mass(deltaT_sec)/mass2;

  double chi = chi_of_spin(spin);
  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i, _traj_.orbital_frequency_isco(chi));
  double t_i = _traj_.time(chi, alpha_i);
  if(duration > -t_i){
    duration = -t_i;
  }

  int time_steps = duration/deltaT + 1;
  return time_steps;
}

class TrajectoryWrapper{
  public:
    TrajectoryWrapper(Vector t, Vector alpha, Vector phase): time(t), alphaOfT(t, alpha), phaseOfT(t, phase) {}
    ~TrajectoryWrapper() {}
    Vector time;
    Spline alphaOfT;
    Spline phaseOfT;
};

TrajectoryWrapper generate_phase_spline(double chi, double phase_i, double alpha_i, double alpha_f){
  int t_steps = 100;
  int power_index = 3;
  double delta_alpha = (pow(alpha_f, power_index) - pow(alpha_i, power_index))/(t_steps - 1);
  Vector t(t_steps);
  Vector phaseOfT(t_steps);
  Vector alphaOfT(t_steps);
  double alpha = pow(alpha_i, power_index);
  for(int j = 0; j < t_steps; j++){
    alphaOfT[j] = pow(alpha, 1./power_index);
    t[j] = _traj_.time(chi, alphaOfT[j]);
    phaseOfT[j] = (_traj_.phase(chi, alphaOfT[j]) - phase_i);
    alpha += delta_alpha;
  }
  return TrajectoryWrapper(t, alphaOfT, phaseOfT);
}

void output_downsampled_trajectory_subfunc(double *t, double *alphaOfT, double *phaseOfT, int t_steps, double chi, double phase_i, double alpha_i, double alpha_f){
  double delta_alpha = (alpha_f - alpha_i)/(t_steps - 1);
  double alpha = alpha_i;
  for(int j = 0; j < t_steps; j++){
    alphaOfT[j] = alpha;
    t[j] = _traj_.time(chi, alpha);
    phaseOfT[j] = (_traj_.phase(chi, alpha) - phase_i);
    alpha += delta_alpha;
  }
}

void output_downsampled_trajectory(double *t, double *alphaOfT, double *phaseOfT, int t_steps, double theta, double phi, double mass1, double mass2, double spin, double r0, double duration_yr, double deltaT_sec){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaT = massratio*seconds_to_solar_mass(deltaT_sec)/mass2;

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  double chi = chi_of_spin(spin);
  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i, _traj_.orbital_frequency_isco(chi));
  double phase_i = _traj_.phase(chi, alpha_i);
  double t_i = _traj_.time(chi, alpha_i);
  if(duration > -t_i){
    duration = -t_i;
  }
  double t_f = t_i + duration;

  double alpha_f = _traj_.orbital_alpha(chi, t_f);
  if(alpha_f < ALPHA_MIN){
    alpha_f = ALPHA_MIN;
  }

  output_downsampled_trajectory_subfunc(t, alphaOfT, phaseOfT, t_steps, chi, phase_i, alpha_i, alpha_f);
}

GWStrain waveform_td(double theta, double phi, double mass1, double mass2, double spin, double r0, double duration_yr, double deltaT_sec, int num_threads){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaT = massratio*seconds_to_solar_mass(deltaT_sec)/mass2;

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  double chi = chi_of_spin(spin);
  double omega_i = kerr_geo_azimuthal_frequency_circ_time(spin, r0);
  double alpha_i = alpha_of_a_omega(spin, omega_i, _traj_.orbital_frequency_isco(chi));
  double phase_i = _traj_.phase(chi, alpha_i);
  double t_i = _traj_.time(chi, alpha_i);
  if(duration > -t_i){
    duration = -t_i;
  }

  double alpha_f = _traj_.orbital_alpha(chi, t_i + duration);
  if(alpha_f < ALPHA_MIN){
    alpha_f = ALPHA_MIN;
  }
  int time_steps = duration/deltaT + 1;
  Vector alpha(time_steps);
  Vector phase(time_steps);
  TrajectoryWrapper trajWrap = generate_phase_spline(chi, phase_i, alpha_i, alpha_f);
  
  alpha[0] = alpha_i;
  phase[0] = 0.;
  Vector alphaComp(time_steps);
  Vector phaseComp(time_steps);
  if(num_threads > 0){
    omp_set_num_threads(num_threads);
  }
  #pragma omp parallel shared(trajWrap, chi, alpha, phase, time_steps, massratio, deltaT, phase_i)
  // #pragma omp parallel shared(_traj_, chi, alpha, phase, time_steps, massratio, deltaT, phase_i)
  {
    #pragma omp for
      for(int j = 1; j < time_steps - 1; j++){
        // alpha[j] = _traj_.orbital_alpha(chi, t_i + deltaT*j);
        // phase[j] = (_traj_.phase(chi, alpha[j]) - phase_i)/massratio;
        // alphaComp[j] = _traj_.orbital_alpha(chi, t_i + deltaT*j);
        // phaseComp[j] = (_traj_.phase(chi, alphaComp[j]) - phase_i)/massratio;
        alpha[j] = trajWrap.alphaOfT.evaluate(t_i + deltaT*j);
        phase[j] = trajWrap.phaseOfT.evaluate(t_i + deltaT*j)/massratio;
      } 
  }
  // for(int j = 1; j < 10; j++){
  //   // alpha[j] = _traj_.orbital_alpha(chi, t_i + deltaT*j);
  //   // phase[j] = (_traj_.phase(chi, alpha[j]) - phase_i)/massratio;
  //   std::cout << "Point i = " << j << "\n";
  //   std::cout << "Alpha difference = " << alphaComp[j] - alpha[j] << "\n";
  //   std::cout << "PhaseComp = " << phaseComp[j] << " and phase = " << phase[j] << "\n";
  //   std::cout << "Phase difference = " << phaseComp[j] - phase[j] << "\n";
  // } 

  // for(int j = time_steps-10; j < time_steps; j++){
  //   // alpha[j] = _traj_.orbital_alpha(chi, t_i + deltaT*j);
  //   // phase[j] = (_traj_.phase(chi, alpha[j]) - phase_i)/massratio;
  //   std::cout << "Point i = " << j << "\n";
  //   std::cout << "Alpha difference = " << alphaComp[j] - alpha[j] << "\n";
  //   std::cout << "PhaseComp = " << phaseComp[j] << " and phase = " << phase[j] << "\n";
  //   std::cout << "Phase difference = " << phaseComp[j] - phase[j] << "\n";
  // } 

  // evaluate at final time of spline to avoid out of bounds error
  // alpha[time_steps - 1] = _traj_.orbital_alpha(chi, t_i + duration);
  // phase[time_steps - 1] = (_traj_.phase(chi, alpha[time_steps - 1]) - phase_i)/massratio;
  alpha[time_steps - 1] = trajWrap.alphaOfT.evaluate(trajWrap.time.back());
  phase[time_steps - 1] = trajWrap.phaseOfT.evaluate(trajWrap.time.back())/massratio;
  GWStrain h(time_steps);
  // std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "It took me " << time_span.count() << " seconds to calculate the phase.";
  // std::cout << std::endl;

  // t1 = std::chrono::high_resolution_clock::now();
  List maxmodes(15);
  harmonic_selection(maxmodes, chi, alpha_i, alpha.back(), REL_MODE_EPS);
  // std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "It took me " << time_span.count() << " seconds.";
  // std::cout << std::endl;
  // t2 = std::chrono::high_resolution_clock::now();
  // time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "It took me " << time_span.count() << " seconds.";
  // std::cout << std::endl;

  // reduced_waveform_harmonic_td(h.plus, h.cross, chi, alpha, maxmodes, theta, phi, phase, num_threads);
  int lmaxmode = maxmodes.size();
  int i = 0;
  while(maxmodes[i] > 0 && i < lmaxmode){
    i++;
  }
  int mmax = i;

  for(int m = 1; m <= mmax ; m++){
    int lmin = (m < 2) ? 2 : m;
    for(int l = lmin; l <= maxmodes[m-1]; l++){
      // t1 = std::chrono::high_resolution_clock::now();
      reduced_waveform_harmonic_td(h.plus, h.cross, chi, alpha, l, m, theta, phi, phase, num_threads);
      // t2 = std::chrono::high_resolution_clock::now();
      // time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      // std::cout << "It took me " << time_span.count() << " seconds.";
      // std::cout << std::endl;
    }
  }

  return h;
}

void generate_waveform_harmonic_td(int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec_temp, int num_threads){
  double deltaT_sec = deltaT_sec_temp;
  if(deltaT_sec == 0.){
    deltaT_sec = 0.9*time_spacing(mass1, mass2, spin, r0, duration_yr);
  }

  std::cout << "Generating waveform with parameters:\n m = "<<mass1<<" M_odot\n M = "<<mass2<<" M_odot\n r0 = "<<r0<<" M\n T = "<<duration_yr<<" yrs\n dt = "<<deltaT_sec<<" sec\n D = "<<dist_Gpc<<" Gpc\n theta = "<<theta<<"\n phi = "<<phi<<"\n";
  GWStrain hTot = waveform_harmonic_td(l, m, theta, phi, mass1, mass2, spin, r0, duration_yr, deltaT_sec, num_threads);
  char buff[100];
  if(deltaT_sec_temp == 0.){
    sprintf(buff, "waveform/h_%d_%d_logm%.4f_logM%.4f_a%.4f_r%.4f_T%.3f_D%.2f_th%.4f_ph%.4f.txt", l, m, log10(mass1), log10(mass2), spin, r0, duration_yr, dist_Gpc, theta, phi);
  }else{
    sprintf(buff, "waveform/h_%d_%d_logm%.4f_logM%.4f_a%.4f_r%.4f_T%.3f_D%.2f_th%.4f_ph%.4f_dt%.1f.txt", l, m, log10(mass1), log10(mass2), spin, r0, duration_yr, dist_Gpc, theta, phi, deltaT_sec);
  }
  std::string filepath = buff;

  std::cout << "Exporting waveform to " << filepath << " \n";
  double rescaledh = scale_strain_amplitude(mass1, dist_Gpc);

  std::ofstream file;
  file.open(filepath);
  file << "t\th+\thx\n";
  for(size_t i = 0; i < hTot.plus.size(); i++){
    sprintf(buff, "%.8e\t%.8e\n", rescaledh*hTot.plus[i], rescaledh*hTot.cross[i]);
    file << buff;
  }

  file.close();
}

void generate_waveform_td(double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec_temp, int num_threads){
  double deltaT_sec = deltaT_sec_temp;
  if(mass1 > mass2){
    double massTemp = mass1;
    mass1 = mass2;
    mass2 = massTemp;
  }
  if(deltaT_sec == 0.){
    deltaT_sec = 0.9*time_spacing(mass1, mass2, spin, r0, duration_yr); // conservative guess on the assumption of max m = 15
  }

  std::cout << "Generating waveform with parameters:\n m = "<<mass1<<" M_odot\n M = "<<mass2<<" M_odot\n MBH spin = "<< spin <<" M\n r0 = "<<r0<<" M\n T = "<<duration_yr<<" yrs\n dt = "<<deltaT_sec<<" sec\n D = "<<dist_Gpc<<" Gpc\n theta = "<<theta<<"\n phi = "<<phi<<"\n";
  GWStrain hTot = waveform_td(theta, phi, mass1, mass2, spin, r0, duration_yr, deltaT_sec, num_threads);
  char buff[100];
  if(deltaT_sec_temp == 0.){
    sprintf(buff, "waveform/h_logm%.4f_logM%.4f_a%.4f_r%.4f_T%.3f_D%.2f_th%.4f_ph%.4f.txt", log10(mass1), log10(mass2), spin, r0, duration_yr, dist_Gpc, theta, phi);
  }else{
    sprintf(buff, "waveform/h_logm%.4f_logM%.4f_a%.4f_r%.4f_T%.3f_D%.2f_th%.4f_ph%.4f_dt%.1f.txt", log10(mass1), log10(mass2), spin, r0, duration_yr, dist_Gpc, theta, phi, deltaT_sec);
  }
  std::string filepath = buff;

  std::cout << "Exporting waveform to " << filepath << " \n";
  double rescaledh = scale_strain_amplitude(mass1, dist_Gpc);

  std::ofstream file;
  file.open(filepath);
  file << "h+\thx\n";
  int hTotSize = hTot.plus.size();
  for(int i = 0; i < hTotSize; i++){
    sprintf(buff, "%.8e\t%.8e\n", rescaledh*hTot.plus[i], rescaledh*hTot.cross[i]);
    file << buff;
  }

  file.close();
}

void output_waveform_td(float *hp, float *hc, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec_temp, int num_threads){
	double deltaT_sec = deltaT_sec_temp;
    if(mass1 > mass2){
      double massTemp = mass1;
      mass1 = mass2;
      mass2 = massTemp;
    }
	if(deltaT_sec == 0.){
		deltaT_sec = 0.9*time_spacing(mass1, mass2, spin, r0, duration_yr);
	}

	GWStrain hTot = waveform_td(theta, phi, mass1, mass2, spin, r0, duration_yr, deltaT_sec, num_threads); // conservative guess on the assumption of max m = 15
    double rescaledh = scale_strain_amplitude(mass1, dist_Gpc);
    int hTotSize = hTot.plus.size();

	for(int i = 0; i < hTotSize; i++){
		hp[i] = rescaledh*hTot.plus[i];
		hc[i] = rescaledh*hTot.cross[i];
	}
}

void output_waveform_harmonic_td(float *hp, float *hc, int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaT_sec_temp, int num_threads){
	double deltaT_sec = deltaT_sec_temp;
    if(mass1 > mass2){
      double massTemp = mass1;
      mass1 = mass2;
      mass2 = massTemp;
    }
	if(deltaT_sec == 0.){
		deltaT_sec = 0.9*time_spacing(mass1, mass2, spin, r0, duration_yr);
	}

	GWStrain hTot = waveform_harmonic_td(l, m, theta, phi, mass1, mass2, spin, r0, duration_yr, deltaT_sec, num_threads); // conservative guess on the assumption of max m = 15
    double rescaledh = scale_strain_amplitude(mass1, dist_Gpc);
    int hTotSize = hTot.plus.size();

	for(int i = 0; i < hTotSize; i++){
		hp[i] = rescaledh*hTot.plus[i];
		hc[i] = rescaledh*hTot.cross[i];
	}
}

////////////////////////////////
// Frequency Domain Waveforms //
////////////////////////////////

void reduced_waveform_harmonic_fd(FloatVector &hp, FloatVector &hc, double chi, const Vector &freq, const int &l, const int &m, const double &t0, const double &theta, const double &phi, const double &phase0, const double &massratio, const double &omegaMin, const double &omegaMax){
  double sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
  double spin = spin_of_chi(chi);
  double domega_df = 2.*M_PI/m;
  double omega = 0., alpha = 0., amp = 0., Phi = 0.;

  double mphi_mod_2pi = fmod(m*phi, 2.*M_PI);
  for(size_t i = 0; i < hp.size(); i++){
    omega = domega_df*freq[i];
    if(omegaMin <= omega && omega <= omegaMax){
      alpha = alpha_of_a_omega(spin, omega, _traj_.orbital_frequency_isco(chi));
      amp = mode_lm(l, m).amplitude(chi, alpha)*sqrt(domega_df*abs(_traj_.time_of_a_omega_derivative(spin, omega))/massratio);
      Phi = mode_lm(l, m).phase(chi, alpha) - fmod(m*(_traj_.phase(chi, alpha) - phase0)/massratio, 2.*M_PI) + mphi_mod_2pi;
      Phi += fmod(m*omega*(_traj_.time(chi, alpha) - t0)/massratio, 2.*M_PI) + 0.25*M_PI;
      hp[i] += amp*sYlm*cos(Phi);
      hc[i] += -amp*sYlm*sin(Phi);
    }
  }
}

void reduced_waveform_harmonic_fd(FloatVector &hpR, FloatVector &hpI, FloatVector &hcR, FloatVector &hcI, double chi, const Vector &freq, const int &l, const int &m, const double &t0, const double &theta, const double &phi, const double &phase0, const double &massratio, const double &omegaMin, const double &omegaMax){
  double sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
  double sYlmMinus = spin_weighted_spherical_harmonic(2, l, m, theta);
  double spin = spin_of_chi(chi);
  double domega_df = 2.*M_PI/m;
  double omega = 0., alpha = 0., amp = 0., Phi = 0.;

  double mphi_mod_2pi = fmod(m*phi, 2.*M_PI);
  for(size_t i = 0; i < hpR.size(); i++){
    omega = domega_df*freq[i];
    if(omegaMin <= omega && omega <= omegaMax){
      alpha = alpha_of_a_omega(spin, omega, _traj_.orbital_frequency_isco(chi));
      amp = mode_lm(l, m).amplitude(chi, alpha)*sqrt(domega_df*abs(_traj_.time_of_a_omega_derivative(spin, omega))/massratio);
      Phi = mode_lm(l, m).phase(chi, alpha) - fmod(m*(_traj_.phase(chi, alpha) - phase0)/massratio, 2.*M_PI) + mphi_mod_2pi;
      Phi += fmod(m*omega*(_traj_.time(chi, alpha) - t0)/massratio, 2.*M_PI) - 0.25*M_PI;
      hpR[i] += 0.5*amp*(sYlm + pow(-1, l + m)*sYlmMinus)*cos(Phi);
      hpI[i] += 0.5*amp*(sYlm + pow(-1, l + m)*sYlmMinus)*sin(Phi);
      hcR[i] += -0.5*amp*(sYlm - pow(-1, l + m)*sYlmMinus)*sin(Phi);
      hcI[i] += 0.5*amp*(sYlm - pow(-1, l + m)*sYlmMinus)*cos(Phi);
    }
  }
}

GWStrainFourier waveform_harmonic_fd(int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double &deltaF_nHz, double theta, double phi){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaF = deltaF_nHz*pow(10, -9)*solar_mass_to_seconds(mass2);

  double chi = chi_of_spin(spin);
  double oISCO = abs(kerr_isco_frequency(spin));

  double omega_i = abs(kerr_geo_azimuthal_frequency_circ_time(spin, r0));
  double alpha_i = alpha_of_a_omega(spin, omega_i, oISCO);
  double t_i = _traj_.time(chi, alpha_i);
  // double phi_i = traj.phase(traj.orbital_alpha(t_i));
  double phi_i = _traj_.phase(chi, alpha_i);

  double t_f = t_i + duration;
  if(t_f > 0){
    t_f = 0;
  }
  double alpha_f = _traj_.orbital_alpha(chi, t_f);
  double omega_f = abs(omega_of_a_alpha(spin, alpha_f, oISCO));

  double freq_i = m*omega_i/(2.*M_PI);
  double freq_f = m*omega_f/(2.*M_PI);

  int freq_steps_skip = freq_i/deltaF;
  while(freq_steps_skip*deltaF > freq_i){
    freq_steps_skip -= 1;
  }
  freq_i = freq_steps_skip*deltaF;
  int freq_steps = 1;
  while(freq_i + (freq_steps - 1)*deltaF < freq_f){
    freq_steps += 1;
  }
  freq_f = freq_i + (freq_steps - 1)*deltaF;

  if(freq_steps < 4){
    std::cout << "(ERROR): Increase frequency resolution to properly sample frequency-domain waveform.\n";
  }

  Vector freq(freq_steps);
  for(int i = 0; i < freq_steps; i++){
    freq[i] = freq_i + i*deltaF;
  }

  // GWStrain h(freq_steps);
  // reduced_waveform_harmonic_fd(h.plus, h.cross, chi, freq, l, m, t_i, theta, phi, phi_i, massratio, omega_i, omega_f);
  GWStrainFourier h(freq_steps);
  reduced_waveform_harmonic_fd(h.plusR, h.plusI, h.crossR, h.crossI, chi, freq, l, m, t_i, theta, phi, phi_i, massratio, omega_i, omega_f);

  return h;
}


GWStrainFourier waveform_fd(double mass1, double mass2, double spin, double r0, double duration_yr, double deltaF_nHz, double theta, double phi){
  double massratio = mass1/mass2;
  double duration = massratio*seconds_to_solar_mass(years_to_seconds(duration_yr))/mass2;
  double deltaF = deltaF_nHz*pow(10, -9)*solar_mass_to_seconds(mass2);

  double chi = chi_of_spin(spin);
  double oISCO = abs(kerr_isco_frequency(spin));

  if(r0 > _traj_.max_orbital_radius(spin)){
    std::cout << "(ERROR): Code only supports max initial distance of r0 = " << _traj_.max_orbital_radius(spin) << " for spin = " << spin << "\n";
  }

  double omega_i = abs(kerr_geo_azimuthal_frequency_circ_time(spin, r0));
  double alpha_i = alpha_of_a_omega(spin, omega_i, oISCO);
  double t_i = _traj_.time(chi, alpha_i);
  // double phi_i = traj.phase(traj.orbital_alpha(t_i));
  double phi_i = _traj_.phase(chi, alpha_i);

  double t_f = t_i + duration;
  if(t_f > 0){
    t_f = 0;
  }
  double alpha_f = _traj_.orbital_alpha(chi, t_f);
  double omega_f = abs(omega_of_a_alpha(spin, alpha_f, oISCO));

    int lmaxmode = 15;
  List maxmodes(lmaxmode);
  harmonic_selection(maxmodes, chi, alpha_i, alpha_f, REL_MODE_EPS);
  int i = 0;
  while(maxmodes[i] > 0 && i < lmaxmode){
    i++;
  }
  int max_m = i;
  double freq_i = omega_i/(2.*M_PI);
  double freq_f = max_m*omega_f/(2.*M_PI);

  int freq_steps_skip = freq_i/deltaF;
  while(freq_steps_skip*deltaF > freq_i){
    freq_steps_skip -= 1;
  }
  freq_i = freq_steps_skip*deltaF;
  int freq_steps = 1;
  while(freq_i + (freq_steps - 1)*deltaF < freq_f){
    freq_steps += 1;
  }
  freq_f = freq_i + (freq_steps - 1)*deltaF;

  if(freq_steps < 4){
    std::cout << "(ERROR): Increase frequency resolution to properly sample frequency-domain waveform.\n";
  }

  Vector freq(freq_steps);
  for(int i = 0; i < freq_steps; i++){
    freq[i] = freq_i + i*deltaF;
  }

  GWStrainFourier h(freq_steps);
  for(int m = 1; m <= max_m; m++){
    int lmin = m < 2 ? 2 : m;
    for(int l = lmin; l <= maxmodes[m-1]; l++){
      reduced_waveform_harmonic_fd(h.plusR, h.plusI, h.crossR, h.crossI, chi, freq, l, m, t_i, theta, phi, phi_i, massratio, omega_i, omega_f);
    }  
  }

  return h;
}

void generate_waveform_harmonic_fd(int l, int m, double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaF_nHz_temp){
  double deltaF_nHz = deltaF_nHz_temp;
  if(deltaF_nHz == 0.){ // based on one over the period
    deltaF_nHz = frequency_spacing(mass1, mass2, spin, r0, duration_yr);
  }
  std::cout << "Generating waveform with parameters:\n m = "<<mass1<<" M_odot\n M = "<<mass2<<" M_odot\n MBH spin = "<<spin<<" M\n r0 = "<<r0<<" M\n T = "<<duration_yr<<" yrs\n df = "<<deltaF_nHz<<" nHz\n D = "<<dist_Gpc<<" Gpc\n theta = "<<theta<<"\n phi = "<<phi<<"\n";
  char buff[100];
  if(deltaF_nHz_temp == 0.){
    sprintf(buff, "waveform/hfd_%d_%d_logm%.4f_logM%.4f_a%.4f_r%.4f_T%.3f_D%.2f_th%.4f_ph%.4f.txt", l, m, log10(mass1), log10(mass2), spin, r0, duration_yr, dist_Gpc, theta, phi);
  }else{
    sprintf(buff, "waveform/hfd_%d_%d_logm%.4f_logM%.4f_a%.4f_r%.4f_T%.3f_D%.2f_th%.4f_ph%.4f_df%.1f.txt", l, m, log10(mass1), log10(mass2), spin, r0, duration_yr, dist_Gpc, theta, phi, deltaF_nHz);
  }
  std::string filepath = buff;
  GWStrainFourier hTot = waveform_harmonic_fd(l, m, mass1, mass2, spin, r0, duration_yr, deltaF_nHz, theta, phi);

  std::cout << "Exporting waveform to " << filepath << " \n";
  double rescaledh = scale_fourier_amplitude(mass1, mass2, dist_Gpc);

  double freq_i = abs(m*kerr_geo_azimuthal_frequency_circ_time(spin, r0)/solar_mass_to_seconds(mass2)/(2.*M_PI));
  int freq_steps_skip = freq_i/deltaF_nHz*pow(10, 9);
  while(freq_steps_skip*deltaF_nHz*pow(10, -9) > freq_i){
    freq_steps_skip -= 1;
  }
  freq_i = freq_steps_skip*deltaF_nHz*pow(10, -9);

  std::ofstream file;
  file.open(filepath);
  file << "f\tRe[h+(f)]\tIm[h+(f)]\tRe[hx(f)]\tIm[hx(f)]\n";
  for(size_t i = 0; i < hTot.plusR.size(); i++){
    sprintf(buff, "%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n", freq_i + i*deltaF_nHz*pow(10, -9), rescaledh*hTot.plusR[i], rescaledh*hTot.plusI[i], rescaledh*hTot.crossR[i], rescaledh*hTot.crossI[i]);
    file << buff;
  }

  file.close();
}

void generate_waveform_fd(double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaF_nHz_temp){

  double deltaF_nHz = deltaF_nHz_temp;
  if(deltaF_nHz == 0.){ // based on one over the period
    deltaF_nHz = frequency_spacing(mass1, mass2, spin, r0, duration_yr);
  }
  std::cout << "Generating waveform with parameters:\n m = "<<mass1<<" M_odot\n M = "<<mass2<<" M_odot\n MBH spin = "<<spin<<" M\n r0 = "<<r0<<" M\n T = "<<duration_yr<<" yrs\n df = "<<deltaF_nHz<<" nHz\n D = "<<dist_Gpc<<" Gpc\n theta = "<<theta<<"\n phi = "<<phi<<"\n";
  GWStrainFourier hTot = waveform_fd(mass1, mass2, spin, r0, duration_yr, deltaF_nHz, theta, phi);
  char buff[100];
  sprintf(buff, "waveform/hfd_logm%.4f_logM%.4f_a%.4f_r%.4f_T%.3f_D%.2f_th%.4f_ph%.4f.txt", log10(mass1), log10(mass2), spin, r0, duration_yr, dist_Gpc, theta, phi);
  std::string filepath = buff;

  // std::string filepath = "waveform/h_logm" + std::to_string(log10(mass1)) + "_logM" + std::to_string(log10(mass2)) + "_r" + std::to_string(r0) + "_T" + std::to_string(duration_yr) + "_dt" + std::to_string(deltaT_sec) + ".txt";
  std::cout << "Exporting waveform to " << filepath << " \n";
  double rescaledh = scale_fourier_amplitude(mass1, mass2, dist_Gpc);

  double freq_i = abs(kerr_geo_azimuthal_frequency_circ_time(spin, r0)/solar_mass_to_seconds(mass2)/(2.*M_PI));
  int freq_steps_skip = freq_i/deltaF_nHz*pow(10, 9);
  while(freq_steps_skip*deltaF_nHz*pow(10, -9) > freq_i){
    freq_steps_skip -= 1;
  }
  freq_i = freq_steps_skip*deltaF_nHz*pow(10, -9);

  std::ofstream file;
  file.open(filepath);
  file << "f\tRe[h+(f)]\tIm[h+(f)]\tRe[hx(f)]\tIm[hx(f)]\n";
  for(size_t i = 0; i < hTot.plusR.size(); i++){
    sprintf(buff, "%.14e\t%.14e\t%.14e\t%.14e\t%.14e\n", freq_i + i*deltaF_nHz*pow(10, -9), rescaledh*hTot.plusR[i], rescaledh*hTot.plusI[i], rescaledh*hTot.crossR[i], rescaledh*hTot.crossI[i]);
    file << buff;
  }

  file.close();
}

//////////////////////////////
// Time-Frequency Waveforms //
//////////////////////////////

void reduced_waveform_harmonic_tf(FloatVector &hp, FloatVector &hc, double chi, const Vector &time, const Vector &freq, const int &l, const int &m, const double &t0, const double &theta, const double &phi, const double &phase0, const double &massratio, const double &omegaMin, const double &omegaMax){
  double sYlm = spin_weighted_spherical_harmonic(-2, l, m, theta);
  double spin = spin_of_chi(chi);
  double domega_df = 2.*M_PI/m;
  double omega = 0., alpha = 0., amp = 0., Phi = 0.;

  double mphi_mod_2pi = fmod(m*phi, 2.*M_PI);
  for(size_t i = 0; i < hp.size(); i++){
    omega = domega_df*freq[i];
    if(omegaMin <= omega && omega <= omegaMax){
      alpha = alpha_of_a_omega(spin, omega, _traj_.orbital_frequency_isco(chi));
      amp = mode_lm(l, m).amplitude(chi, alpha);
      Phi = mode_lm(l, m).phase(chi, alpha) - fmod(m*(_traj_.phase(chi, alpha) - phase0)/massratio, 2.*M_PI) + mphi_mod_2pi;
      Phi += fmod(m*omega*(_traj_.time(chi, alpha) - t0)/massratio, 2.*M_PI) - 0.25*M_PI;
      hp[i] += amp*sYlm*cos(Phi);
      hc[i] += -amp*sYlm*sin(Phi);
    }
  }
}
