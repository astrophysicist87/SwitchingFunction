#ifndef THERMO_EXII_H
#define THERMO_EXII_H


double Total_Pressure_exII(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList, double epsilon0 = -1.0,
 double Press_exII_guess = -1.0);

double Total_Energy_Dens_exII(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList, double epsilon0 = -1.0,
 double Press_exII_guess = -1.0);

double Total_Baryon_Dens_exII(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList, double epsilon0 = -1.0,
 double Press_exII_guess = -1.0);

double Total_Number_Dens_exII(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList, double epsilon0 = -1.0,
 double Press_exII_guess = -1.0);

double Total_Entropy_Dens_exII(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList, double epsilon0 = -1.0,
 double Press_exII_guess = -1.0);










#endif