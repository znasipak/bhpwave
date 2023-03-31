#include "snr.hpp"

double sensitivity_curve_LISA(double f){
  double f1 = 0.4e-3;
  double f2 = 25.e-3;
  double SI = (5.76e-48)*(1. + pow(f1/f, 2));
  double SII = 3.6e-41;
  double R = 1. + pow(f/f2, 2);

  return 10./3.*(SI*pow(2.*M_PI*f, -4) + SII)*R;
}

Vector sensitivity_curve_LISA(Vector f){
  double f1 = 0.4e-3;
  double f2 = 25.e-3;
  double SIcoeff = 5.76e-48;
  double SIIcoeff = 3.6e-41;
  Vector Sn(f.size());

  for(size_t i = 0; i < Sn.size(); i++){
    double SI = SIcoeff*(1. + pow(f1/f[i], 2));
    double SII = SIIcoeff;
    double R = 1. + pow(f[i]/f2, 2);
    Sn[i] = 10./3.*(SI*pow(2.*M_PI*f[i], -4) + SII)*R;
  }

  return Sn;
}

double snr_LISA(double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaF_nHz_temp){
  double deltaF_nHz = deltaF_nHz_temp;
  if(deltaF_nHz == 0.){ // based on one over the period
    deltaF_nHz = frequency_spacing(mass1, mass2, spin, r0, duration_yr);
  }

  std::cout << "Generating waveform with parameters:\n m = "<<mass1<<" M_odot\n M = "<<mass2<<" M_odot\n MBH spin = "<<spin<<" M\n r0 = "<<r0<<" M\n T = "<<duration_yr<<" yrs\n df = "<<deltaF_nHz<<" nHz\n D = "<<dist_Gpc<<" Gpc\n theta = "<<theta<<"\n phi = "<<phi<<"\n";

  GWStrain hTot = waveform_fd(mass1, mass2, spin, r0, duration_yr, deltaF_nHz, theta, phi);

  double rescaledh = scale_fourier_amplitude(mass1, mass2, dist_Gpc);

  double freq_i = abs(kerr_geo_azimuthal_frequency_circ_time(spin, r0)/solar_mass_to_seconds(mass2)/(2.*M_PI));
  int freq_steps_skip = freq_i/deltaF_nHz*pow(10, 9);
  while(freq_steps_skip*deltaF_nHz*pow(10, -9) > freq_i){
    freq_steps_skip -= 1;
  }
  freq_i = freq_steps_skip*deltaF_nHz*pow(10, -9);

  double snr2 = 0.;
  for(size_t i = 0; i < hTot.plus.size(); i++){
    snr2 += 4.*(pow(hTot.plus[i], 2) + pow(hTot.cross[i], 2))/sensitivity_curve_LISA(freq_i + i*deltaF_nHz*pow(10, -9));
  }

  return rescaledh*sqrt(snr2*deltaF_nHz*pow(10, -9));
}
