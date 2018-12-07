#ifndef THERMO_PT_H
#define THERMO_PT_H

#include <vector>
#include "particle.hpp"


//Fermi/Bose/Classical stats distribution function
double equilib_distfcn(double Energy, double Temp, double mu_particle, 
 QUANTUM_TYPE quantumtype);

//-----------------------------------------------------------------------------
//Code to compute (partial) thermal properties of one species of gas:

double Partial_Pressure_pt(double Temp, double mu_particle, 
 const Particle& particle);

double Partial_Number_Dens_pt(double Temp, double mu_particle,
 const Particle& particle);

double Partial_Baryon_Dens_pt(double Temp, double mu_particle,
 const Particle& particle);

double Partial_Energy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle);

double Partial_Entropy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle);

//-----------------------------------------------------------------------------
//Code to compute thermal properties of a point particle gas of many particle 
//species

double Total_Pressure_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList);

double Total_Energy_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList);

double Total_Number_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList);

double Total_Baryon_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList);

double Total_Entropy_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList);

//-----------------------------------------------------------------------------







#endif
