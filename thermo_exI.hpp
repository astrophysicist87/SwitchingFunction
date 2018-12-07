#ifndef THERMO_EXI_H
#define THERMO_EXI_H

#include <vector>
#include "particle.hpp"

//-----------------------------------------------------------------------------
//Code to switch between Temp, mu_B and T_star, mu_B_star for exI model

double Temp_exI(double T_star, double mu_B_star, double epsilon0, 
 const std::vector<Particle>& ParticleList);

double mu_B_exI(double T_star, double mu_B_star, double epsilon0,
 const std::vector<Particle>& ParticleList);

void Find_Temp_star_mu_B_star_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList, double& T_star, double& mu_B_star);


//-----------------------------------------------------------------------------
//Code to compute thermal properties of a point particle gas of many particle 
//species

double Total_Pressure_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList);

double Total_Energy_Dens_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList);

double Total_Number_Dens_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList);

double Total_Baryon_Dens_exI(double Temp, double mu_B, double epsilon0, 
 const std::vector<Particle>& ParticleList);

double Total_Entropy_Dens_exI(double Temp, double mu_B, double epsilon0, 
 const std::vector<Particle>& ParticleList);

//-----------------------------------------------------------------------------







#endif
