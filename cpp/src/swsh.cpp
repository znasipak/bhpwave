// Spin-weighted spherical harmonics

#include "swsh.hpp"

double factorial(int n){
    return gsl_sf_fact(n);
}

double log_factorial(int n){
    return gsl_sf_lnfact(n);
}

double binomial(int n, int m){
    return gsl_sf_choose(n, m);
}

double spin_weighted_harmonic_prefactor(int s, int l, int m){
    double factorial_factor = exp(log_factorial(l + m) + log_factorial(l - m) - log_factorial(l + s) - log_factorial(l - s));
    return pow(-1, m)*sqrt(factorial_factor*(2*l + 1)/(4.*M_PI));
}

double spherical_harmonic(const int &l, const int &m, const double &th){
	if( m < 0 && l >= abs(m) ){
		return pow(-1, m)*gsl_sf_legendre_sphPlm(l, -m, cos(th));
	}else if(l >= abs(m)){
		return gsl_sf_legendre_sphPlm(l, m, cos(th));
	}else{
		return 0.;
	}
}

Vector spherical_harmonic(const int &l, const int &m, const Vector &th){
	Vector ylm(th.size());
    for(size_t i = 0; i < th.size(); i++){
        ylm[i] = spherical_harmonic(l, m, th[i]);
    }
    return ylm;
}

Complex spherical_harmonic(const int &l, const int &m, const double &th, const double &ph){
	return spherical_harmonic(l, m, th)*exp(Complex(0., m*ph));
}

double spin_weighted_sum(int s, int l, int m, double z){
    int rmax = l - s;
    double yslm = 0.;
    for(int r = 0; r <= rmax; r++){
        if(r + s - m >= 0){
            double term = binomial(l - s, r)*binomial(l + s, r + s - m);
            term *= pow(z - 1., rmax - r)*pow(z + 1., r);
            yslm += term;
        }
    }
    return yslm;
}

// double spin_weighted_sum(int s, int l, int m, double z){
//     int rmax = l - s;
//     double yslm = 0.;
//     int r = 0;
//     double term = binomial(l + s, s - m);
//     term *= pow(z - 1., rmax);
//     yslm += term;
//     for(r = 1; r < rmax; r++){
//         if(r + s - m >= 0){
//             term = binomial(l - s, r)*binomial(l + s, r + s - m);
//             term *= pow(z - 1., rmax - r)*pow(z + 1., r);
//             yslm += term;
//         }
//     }
//     r = rmax;
//     if(r + s - m >= 0){
//         term = binomial(l - s, r)*binomial(l + s, r + s - m);
//         term *= pow(z + 1., r);
//         yslm += term;
//     }

//     return yslm;
// }

double spin_weighted_sum_dz(int s, int l, int m, double z){
    int rmax = l - s;
    double yslm = 0.;
    int r = 0;
    double term = binomial(l + s, s - m);
    term *= rmax*pow(z - 1., rmax - 1);
    yslm += term;
    for(r = 1; r < rmax; r++){
        if(r + s - m >= 0){
            term = binomial(l - s, r)*binomial(l + s, r + s - m);
            term *= (rmax - r)*pow(z - 1., rmax - r - 1)*pow(z + 1., r) + r*pow(z - 1., rmax - r)*pow(z + 1., r - 1);
            yslm += term;
        }
    }
    r = rmax;
    if(r + s - m >= 0){
        term = binomial(l - s, r)*binomial(l + s, r + s - m);
        term *= r*pow(z + 1., r - 1);
        yslm += term;
    }

    return yslm;
}

double spin_weighted_spherical_harmonic(int s, int l, int m, double theta){
  if(s == 0){
    return spherical_harmonic(l, m, theta);
  }else if(s + m < 0){
    return pow(-1, s+m)*spin_weighted_spherical_harmonic(-s, l, -m, theta);
  }else if(theta > 0.5*M_PI){
    return pow(-1, l + m)*spin_weighted_spherical_harmonic(-s, l, m, M_PI - theta);
  }
  double z = cos(theta);
  double pref = pow(0.5, l)*spin_weighted_harmonic_prefactor(s, l, m)*pow(1. - z, 0.5*(s + m))*pow(1. + z, 0.5*(s - m));
  return pref*spin_weighted_sum(s, l, m, z);
}

double spin_weighted_spherical_harmonic_dz(int s, int l, int m, double theta){
  if(s + m < 0){
    return pow(-1, s+m)*spin_weighted_spherical_harmonic_dz(-s, l, -m, theta);
  }else if(theta > 0.5*M_PI){
    return -pow(-1, l + m)*spin_weighted_spherical_harmonic_dz(-s, l, m, M_PI - theta);
  }
  double z = cos(theta);
  double pref = pow(0.5, l)*spin_weighted_harmonic_prefactor(s, l, m)*pow(1. - z, 0.5*(s + m))*pow(1. + z, 0.5*(s - m));
  double dpref = pow(0.5, l)*spin_weighted_harmonic_prefactor(s, l, m);
  dpref *= -0.5*(s + m)*pow(1. - z, 0.5*(s + m) - 1)*pow(1. + z, 0.5*(s - m)) + 0.5*(s - m)*pow(1. - z, 0.5*(s + m))*pow(1. + z, 0.5*(s - m) - 1);
  return dpref*spin_weighted_sum(s, l, m, z) + dpref*spin_weighted_sum_dz(s, l, m, z);
}

Vector spin_weighted_spherical_harmonic(int st, int l, int mt, Vector theta){
  if(st == 0){
    return spherical_harmonic(l, mt, theta);
  }

  double pref;
  int s = st;
  int m = mt;
  if(s + m < 0){
    s = -st;
    m = -mt;
    pref = pow(0.5, l)*pow(-1, s+m)*spin_weighted_harmonic_prefactor(s, l, m);
  }else{
    pref = pow(0.5, l)*spin_weighted_harmonic_prefactor(s, l, m);
  }

  Vector yslm(theta.size());
  double z;
  for(size_t i = 0; i < theta.size(); i++){
    if(theta[i] > 0.5*M_PI){
      z = -cos(theta[i]);
      yslm[i] = pow(-1, l + m)*pow(1. - z, 0.5*(s + m))*pow(1. + z, 0.5*(s - m))*pref*spin_weighted_sum(-s, l, m, z);
    }else{
      z = cos(theta[i]);
      yslm[i] = pow(1. - z, 0.5*(s + m))*pow(1. + z, 0.5*(s - m))*pref*spin_weighted_sum(s, l, m, z);
    }
  }

  return yslm;
}

void spin_weighted_spherical_harmonic(double *yslm, int pts_num, int st, int l, int mt, double *theta){
  if(st == 0){
    for(int i = 0; i < pts_num; i++){
      yslm[i] = spherical_harmonic(l, mt, theta[i]);
    }
  }

  double pref;
  int s = st;
  int m = mt;
  double z;
  
  if(s + m < 0){
    s = -st;
    m = -mt;
    pref = pow(0.5, l)*pow(-1, s+m)*spin_weighted_harmonic_prefactor(s, l, m);
  }else{
    pref = pow(0.5, l)*spin_weighted_harmonic_prefactor(s, l, m);
  }

  for(size_t i = 0; i < pts_num; i++){
    if(theta[i] > 0.5*M_PI){
      z = -cos(theta[i]);
      yslm[i] = pow(-1, l + m)*pow(1. - z, 0.5*(s + m))*pow(1. + z, 0.5*(s - m))*pref*spin_weighted_sum(-s, l, m, z);
    }else{
      z = cos(theta[i]);
      yslm[i] = pow(1. - z, 0.5*(s + m))*pow(1. + z, 0.5*(s - m))*pref*spin_weighted_sum(s, l, m, z);
    }
  }
}

Complex spin_weighted_spherical_harmonic(int s, int l, int m, double theta, double phi){
  return spin_weighted_spherical_harmonic(s, l, m, theta)*exp(Complex(0., m*phi));
}