#include "snr.hpp"

void save_data(std::string filename, std::string header, Matrix data){
  std::ofstream file;
	file.open(filename);

  char buff[20];
  file << header << "\n";
  for(size_t i = 0; i < data[0].size(); i++){
    for(size_t j = 0; j < data.size(); j++){
      if(j < data.size() - 1){
        sprintf(buff, "%.15e\t", data[j][i]);
      }else{
        sprintf(buff, "%.15e\n", data[j][i]);
      }
      file << buff;
    }
  }

	file.close();
}

int main(){
  // TrajectoryStruct trajStruct = read_trajectory_data();
  // for(size_t i = 0; i < trajStruct.alpha.size(); i++){
  //   std::cout << trajStruct.alpha[i] << "\n";
  // }
  // for(size_t i = 0; i < trajStruct.chi.size(); i++){
  //   std::cout << trajStruct.chi[i] << "\n";
  // }

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  load_waveform_data();
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "Took " << time_span.count() << " seconds to load data. \n";

  // TrajectorySpline2D traj(read_trajectory_data());
  // double omega0 = kerr_geo_azimuthal_frequency_circ_time(0.9, 20.);
  // std::cout <<  traj.time_of_a_omega(0.9, omega0) << "\n";

  // TrajectoryStruct trajStruct = read_trajectory_data();
  // for(size_t i = 0; i < trajStruct.tMax.size(); i++){
  //   std::cout << trajStruct.tMax[i] << "\n";
  // }
  // for(size_t i = 0; i < trajStruct.tMax.size(); i++){
  //   std::cout << trajStruct.chi[i] << "\n";
  // }
  // TrajectorySpline2D traj(read_trajectory_data());
  // double omega0 = kerr_geo_azimuthal_frequency_circ_time(0., 10.);
  // // double omega0 = 0.031682894767489804;
  // double alpha0 = alpha_of_a_omega(0., omega0);
  // double chi0 = chi_of_spin(0.);
  // // std::cout << omega0 << "\n";
  // // std::cout << alpha0 << "\n";
  // double t0 = traj.time_of_a_omega(0., omega0);
  // double t0_2 = traj.time(chi0, alpha0);
  // std::cout << "Min orbital frequency = " << traj.min_orbital_frequency(0.) << "\n";
  // std::cout << "Max orbital frequency = " << traj.max_orbital_frequency(0.) << "\n";
  // std::cout << "Time comparison: t0 = " << abs(1. - t0/t0_2) << "\n";
  // std::cout << "Freq comparison: omega0 = " << abs(1. - omega0/traj.orbital_frequency(0., t0)) << "\n";
  // std::cout << "Alpha comparison: alpha0 = " << abs(1. - alpha0/traj.orbital_alpha(chi0, t0)) << "\n";
  // std::cout << "Max comparison: max_time_before_merger = " << abs(1. - traj.max_time_before_merger(0.)/traj.time_of_a_omega(0., traj.min_orbital_frequency(0.))) << "\n";
  // std::cout << "Min comparison: min_orbital_frequency = " << abs(1. - traj.min_orbital_frequency(0.)/traj.orbital_frequency(0., traj.max_time_before_merger(0.))) << "\n";
  // std::cout << "Min comparison: min_alpha_frequency = " << abs(1. - traj.orbital_alpha(chi0, traj.max_time_before_merger(0.))) << "\n";
  //
  // double Phi0 = traj.phase(1., traj.orbital_alpha(1., t0));
  // int timeSampleNumber = 100001;

  // double dt = 90.97841171825908/1.e5;
  // while(t0 + dt*(timeSampleNumber - 1) > 0.){
  //   timeSampleNumber -= 1;
  // }
  // std::cout << 1.e5*(t0 + dt*(timeSampleNumber - 1)) << "\n";
  // Matrix data(4, Vector(timeSampleNumber, 0.));
  // for(int i = 0; i < timeSampleNumber; i++){
  //   data[0][i] = 1.e5*i*dt;
  //   data[1][i] = 2*1.e5*(traj.phase(1., traj.orbital_alpha(1., t0 + i*dt)) - Phi0);
  //   data[2][i] = traj.orbital_frequency(0., t0 + i*dt);
  //   data[3][i] = traj.orbital_frequency_derivative(0., t0 + i*dt);
  // }
  // save_data("./post_processing/phase_hughes_circ_v2.txt", "t\tPhi_m\tOmega\tOmegaDot", data);

  // double dt = 90.80578054657609/1.e5;
  // double dt = 90.7634694189251/1.e5;
  // while(t0 + dt*(timeSampleNumber - 1) > 0.){
  //   timeSampleNumber -= 1;
  // }
  // Matrix data(4, Vector(timeSampleNumber, 0.));
  // for(int i = 0; i < timeSampleNumber; i++){
  //   data[0][i] = 1.e5*i*dt;
  //   data[1][i] = 1.e5*(traj.phase(1., traj.orbital_alpha(1., t0 + i*dt)) - Phi0);
  //   data[2][i] = traj.orbital_frequency(0., t0 + i*dt);
  //   data[3][i] = 1.e-5*traj.orbital_frequency_derivative(0., t0 + i*dt);
  // }
  // save_data("./post_processing/phase_hughes_dense_v2.txt", "t\tPhi_m\tOmega\tOmegaDot", data);
  //
  // dt = 90.80622839774117/1.e5;
  // for(int i = 0; i < timeSampleNumber; i++){
  //   data[0][i] = 1.e5*i*dt;
  //   data[1][i] = 2*1.e5*(traj.phase(1., traj.orbital_alpha(1., t0 + i*dt)) - Phi0);
  //   data[2][i] = traj.orbital_frequency(0., t0 + i*dt);
  //   data[3][i] = traj.orbital_frequency_derivative(0., t0 + i*dt);
  // }
  // save_data("./post_processing/phase_hughes_ecc.txt", "t\tPhi_m\tOmega\tOmegaDot", data);
  //
  // dt = 90.80770966193103/1.e5;
  // for(int i = 0; i < timeSampleNumber; i++){
  //   data[0][i] = 1.e5*i*dt;
  //   data[1][i] = 2*1.e5*(traj.phase(1., traj.orbital_alpha(1., t0 + i*dt)) - Phi0);
  //   data[2][i] = traj.orbital_frequency(0., t0 + i*dt);
  //   data[3][i] = traj.orbital_frequency_derivative(0., t0 + i*dt);
  // }
  // save_data("./post_processing/phase_hughes_ecc2.txt", "t\tPhi_m\tOmega\tOmegaDot", data);

  t1 = std::chrono::high_resolution_clock::now();
  generate_waveform_td(10., 1.e6, 0.1, 10., 10., 2., 0.1, 0.1);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "Took " << time_span.count() << " seconds to generate waveform. \n";

  t1 = std::chrono::high_resolution_clock::now();
  generate_trajectory_td(2., 1.e6, 0.995, 6., 0., 1.5);
  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "Took " << time_span.count() << " seconds to generate trajectory. \n";

  // t1 = std::chrono::high_resolution_clock::now();
  // generate_waveform_fd(2., 1.e6, 0.995, 6., 1.5, 2., 0.5*M_PI, 0.);
  // t2 = std::chrono::high_resolution_clock::now();
  // time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "Took " << time_span.count() << " seconds to generate waveform. \n";

  // t1 = std::chrono::high_resolution_clock::now();
  // generate_waveform_td(2., 1.e6, 0.9, 6., 1.5, 2., 0.3*M_PI, 0.);
  // t2 = std::chrono::high_resolution_clock::now();
  // time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "Took " << time_span.count() << " seconds to generate waveform. \n";

  // t1 = std::chrono::high_resolution_clock::now();
  // generate_waveform_fd(2., 1.e6, 0.9, 6., 1.5, 2., 0.3*M_PI, 0.);
  // t2 = std::chrono::high_resolution_clock::now();
  // time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "Took " << time_span.count() << " seconds to generate waveform. \n";

  // t1 = std::chrono::high_resolution_clock::now();
  // generate_waveform_td(10., 1.e6, 0.9, 6., 1.5, 2., 0.3*M_PI, 0.);
  // t2 = std::chrono::high_resolution_clock::now();
  // time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "Took " << time_span.count() << " seconds to generate waveform. \n";

  // t1 = std::chrono::high_resolution_clock::now();
  // generate_waveform_fd(10., 1.e6, 0.9, 6., 1.5, 2., 0.3*M_PI, 0.);
  // t2 = std::chrono::high_resolution_clock::now();
  // time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  // std::cout << "Took " << time_span.count() << " seconds to generate waveform. \n";

  return 0;
}
