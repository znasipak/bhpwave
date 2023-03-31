#ifndef SNR_HPP
#define SNR_HPP

#include "waveform.hpp"

double sensitivity_curve_LISA(double f);
Vector sensitivity_curve_LISA(Vector f);

double snr_LISA(double mass1, double mass2, double spin, double r0, double duration_yr, double dist_Gpc, double theta, double phi, double deltaF_nHz_temp = 0.);

#endif